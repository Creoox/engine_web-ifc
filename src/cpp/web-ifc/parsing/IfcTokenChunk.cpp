/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
 
#include <array>
#include <spdlog/spdlog.h>
#include "IfcTokenStream.h"

// Character-class dispatch table
// Maps every byte value (0-255) to a dispatch tag in O(1): one array lookup,
// no comparisons, no range tests.  Fits in 256 bytes = 4 cache lines.
namespace {
  enum CharClass : uint8_t {
	CC_OTHER = 0,	// unrecognised byte - silently advance (also the default zero-init)
	CC_WS,			// space  \t  \n  \r
	CC_STR,			// '
	CC_REF,			// #
	CC_EMPTY,		// $
	CC_STAR,		// *
	CC_SETB,		// (
	CC_SETE,		// )
	CC_LINE,		// ;
	CC_DOT,			// .
	CC_DIGIT,		// 0-9
	CC_ALPHA,		// A-Z  a-z
	CC_CTRL,		// control chars < 32, excluding whitespace
  };

  constexpr std::array<uint8_t, 256> makeCharClassTable()
  {
	std::array<uint8_t, 256> t{};			// zero-init → CC_OTHER everywhere

	for (int i = 0; i < 32; ++i) t[i] = CC_CTRL;   // all controls are errors …
	t['\t'] = CC_WS; t['\n'] = CC_WS;			   // … except whitespace
	t['\r'] = CC_WS; t[' ']  = CC_WS;

	t['\''] = CC_STR;   t['#'] = CC_REF;   t['$'] = CC_EMPTY;
	t['*']  = CC_STAR;  t['('] = CC_SETB;  t[')'] = CC_SETE;
	t[';']  = CC_LINE;  t['.'] = CC_DOT;

	for (int c = '0'; c <= '9'; ++c) t[c] = CC_DIGIT;
	for (int c = 'A'; c <= 'Z'; ++c) t[c] = CC_ALPHA;
	for (int c = 'a'; c <= 'z'; ++c) t[c] = CC_ALPHA;

	return t;
  }
  constexpr auto kCharClass = makeCharClassTable();
}

namespace webifc::parsing
{
	IfcTokenStream::IfcTokenChunk::IfcTokenChunk(const size_t chunkSize, const size_t startRef, const size_t fileStartRef, IfcFileStream* fileStream) : _startRef(startRef), _fileStartRef(fileStartRef), _chunkSize(chunkSize), _fileStream(fileStream)
	{
		_chunkData = nullptr;
		_loaded = true;
		_currentSize = 0;
		if (_fileStream != nullptr) Load();
	}

	bool IfcTokenStream::IfcTokenChunk::Clear(bool force)
	{
		if (_fileStream == nullptr && !force) return false;
		if (_chunkData != nullptr) { delete[] _chunkData; _chunkData = nullptr; }
		_loaded = false;
		return true;
	}

	bool IfcTokenStream::IfcTokenChunk::Clear()
	{
		return Clear(false);
	}

	size_t IfcTokenStream::IfcTokenChunk::GetTokenRef()
	{
		return _startRef;
	}

	size_t IfcTokenStream::IfcTokenChunk::TokenSize()
	{
		return _currentSize;
	}

	bool IfcTokenStream::IfcTokenChunk::IsLoaded()
	{
		return _loaded;
	}

	size_t IfcTokenStream::IfcTokenChunk::GetMaxSize()
	{
		return _chunkSize;
	}

	std::string_view IfcTokenStream::IfcTokenChunk::ReadString(const size_t ptr, const size_t size)
	{
		if (!_loaded) Load();
		return std::string_view((char*)_chunkData + ptr, size);
	}

	void IfcTokenStream::IfcTokenChunk::Push(void* v, const size_t size)
	{
		if (_chunkData == nullptr) _chunkData = new uint8_t[_chunkSize];
		_currentSize += size;
		if (_currentSize > _chunkSize) {
			size_t newCapacity = _currentSize * 2;
			uint8_t* tmp = _chunkData;
			_chunkData = new uint8_t[newCapacity];
			std::memcpy(_chunkData, tmp, _currentSize - size);
			_chunkSize = newCapacity;
			delete[] tmp;
		}
		std::memcpy(_chunkData + _currentSize - size, v, size);
	}

	void IfcTokenStream::IfcTokenChunk::Load()
	{
		_chunkData = new uint8_t[_chunkSize];
		_loaded = true;
		if (_fileStream->GetRef() != _fileStartRef) _fileStream->Go(_fileStartRef);
		std::vector<char> temp;
		temp.reserve(50);
		_currentSize = 0;
		while (!_fileStream->IsAtEnd() && _currentSize < _chunkSize)
		{
			const char c = _fileStream->Get();

			switch (kCharClass[(unsigned char)c])
			{
				// Whitespace
			case CC_WS:
				_fileStream->Forward();
				continue;   // back to top of while - no outer Forward() needed

				// String literal
			case CC_STR:
			{
				_fileStream->Forward();
				temp.clear();
				// apparently I dont fully understand strings in IFC yet
				// this example from uptown shows that escaping is not used: 'Type G5 - 800kg/m\X2\00B2\X0\';
				// this example from revit shows that double quotes are used as one quote: 'RPC Tree - Deciduous:Scarlet Oak - 42'':946835'
				// turns out this is just part of ISO 10303-21, thanks ottosson!
				while (true)
				{
					const char sc = _fileStream->Get(); // read once, reuse below
					temp.push_back(sc);
					// if its a quote, maybe its the end of the string
					if (sc == '\'')
					{
						// if there's another quote behind it, its not
						_fileStream->Forward();
						if (_fileStream->Get() == '\'')
						{
							// we also bump pos, otherwise we still break next loop...
							temp.push_back(_fileStream->Get());
						}
						else
						{
							_fileStream->Back();
							temp.pop_back();
							break;
						}
					}
					_fileStream->Forward();
				}
				Push<uint8_t>(IfcTokenType::STRING);
				Push<uint16_t>(temp.size());
				if (temp.size() > 0) Push(temp.data(), temp.size());
				break;
			}

			// Entity reference (#nnn)
			case CC_REF:
			{
				_fileStream->Forward();
				uint32_t num = 0;
				char d = _fileStream->Get();
				while (d >= '0' && d <= '9')
				{
					num = num * 10 + (d - '0');
					_fileStream->Forward();
					d = _fileStream->Get();
				}
				Push<uint8_t>(IfcTokenType::REF);
				Push<uint32_t>(num);
				continue;   // stream already past the last digit
			}

			// Empty / unset ($)
			case CC_EMPTY:
				Push<uint8_t>(IfcTokenType::EMPTY);
				break;

				// Comment (/*…*/) or UNKNOWN (*)
			case CC_STAR:
				if (_fileStream->Prev() == '/')
				{
					_fileStream->Forward();
					while (!(_fileStream->Prev() == '*' && _fileStream->Get() == '/'))
						_fileStream->Forward();
				}
				else Push<uint8_t>(IfcTokenType::UNKNOWN);
				break;

				// Punctuation
			case CC_SETB:  Push<uint8_t>(IfcTokenType::SET_BEGIN);  break;
			case CC_SETE:  Push<uint8_t>(IfcTokenType::SET_END);	break;
			case CC_LINE:  Push<uint8_t>(IfcTokenType::LINE_END);   break;

				// Number (integer or real)
			case CC_DIGIT:
			{
				// TODO: consider storing the number in binary form instead of as a string, to save space and parsing time later.  
				// The main obstacle is that we need to preserve the exact formatting of the number, including leading zeros, decimal point,
				// exponent marker, etc., because these may be significant for some applications.
				// If we convert to binary and back to string later, we may lose this formatting information.
				// One possible solution is to store both the original string and the parsed binary value, but this would increase memory usage.  
				temp.clear();
				if (_fileStream->Prev() == '-') temp.push_back('-');
				char d = _fileStream->Get();
				bool isFrac = false;
				while ((d >= '0' && d <= '9') || d == '.' || d == 'e' || d == 'E' || d == '-' || d == '+')
				{
					temp.push_back(d);
					if (d == '.' || d == 'E') isFrac = true;
					_fileStream->Forward();
					d = _fileStream->Get();
				}
				Push<uint8_t>(isFrac ? IfcTokenType::REAL : IfcTokenType::INTEGER);
				Push<uint16_t>(temp.size());
				Push(temp.data(), temp.size());
				continue;   // stream already past the last digit
			}

			// Enumeration (.NAME.)
			case CC_DOT:
			{
				temp.clear();
				_fileStream->Forward();
				char d = _fileStream->Get();
				while (d != '.')
				{
					temp.push_back(d);
					_fileStream->Forward();
					d = _fileStream->Get();
				}
				Push<uint8_t>(IfcTokenType::ENUM);
				Push<uint16_t>(temp.size());
				Push(temp.data(), temp.size());
				break;
			}

			// Keyword / type label
			case CC_ALPHA:
			{
				temp.clear();
				char d = _fileStream->Get();
				while ((d >= 'A' && d <= 'Z') || (d >= 'a' && d <= 'z') || (d >= '0' && d <= '9') || d == '_')
				{
					temp.push_back(d);
					_fileStream->Forward();
					d = _fileStream->Get();
				}
				Push<uint8_t>(IfcTokenType::LABEL);
				Push<uint16_t>(temp.size());
				Push(temp.data(), temp.size());
				continue;   // stream already past the last label char
			}

			// Invalid control character (e.g. embedded '\0')
			case CC_CTRL:
			{
				// CC_WS already handled above; this branch fires only for genuine control bytes such as '\0'.
				spdlog::error("[IfcTokenStream] Invalid control character 0x{:02X} at file position {}; skipping to next statement",
					static_cast<uint8_t>(c), _fileStream->GetRef());
				// Advance PAST the invalid byte, then scan forward until a safe
				// statement delimiter.  Guard IsAtEnd() after each Forward().
				while (!_fileStream->IsAtEnd())
				{
					_fileStream->Forward();
					if (_fileStream->IsAtEnd()) break;
					const char skip = _fileStream->Get();
					if (skip == ';') { Push<uint8_t>(IfcTokenType::LINE_END); break; }
					if (skip == '\n') break;
				}
				break;
			}

			// CC_OTHER (= 0) has no case: the switch drops through to the outer
			// _fileStream->Forward() below, silently consuming unrecognised bytes
			// (e.g. '-' that precedes a number, '/' in comments, etc.).
			}

			_fileStream->Forward();
		}
	}
}