/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
 
#pragma once
 
#include <vector>
#include <istream>
#include <iostream>
#include <functional>
#include <string_view>
#include <cstring>
#include <cstdint>
 
namespace webifc::parsing
{
  
  enum IfcTokenType : char
  {
    UNKNOWN = 0,
    STRING,
    LABEL,
    ENUM,
    REAL,
    REF,
    EMPTY,
    SET_BEGIN,
    SET_END,
    LINE_END,
    INTEGER
  };
  
  
  class IfcTokenStream 
  {
      public:
        IfcTokenStream(const size_t chunkSize, const uint64_t maxChunks);
        ~IfcTokenStream();
        void SetTokenSource(const std::function<uint32_t(char *, size_t, size_t)> &requestData);
        void SetTokenSource(std::istream &requestData);
        template <typename T> T Read()
        {
          if (!_cChunk->IsLoaded()) {
            checkMemory();
            _activeChunks++;
          }
          T v =  _cChunk->Read<T>(_readPtr);
          Forward(sizeof(T));
          return v;
        }
        template <typename T> void Push(T input)
        {
          Push(&input,sizeof(T));
        }
        void Push(void *v, const size_t size);
        void Forward(const size_t size);
        std::string_view ReadString();
        void Back();
        bool IsAtEnd();
        void MoveTo(const size_t pos);
        size_t GetReadOffset();
        size_t GetTotalSize();
        IfcTokenStream * Clone();
        size_t GetNoLines();

      private:
        void checkMemory();
        size_t _readPtr = 0;
      	size_t _currentChunk = 0;
        size_t _activeChunks = 0;
        size_t _chunkSize;
        uint64_t _maxChunks;
        class IfcFileStream
        {
          public:
            IfcFileStream(const std::function<uint32_t(char *, size_t, size_t)> &requestData, const uint32_t size);
            ~IfcFileStream();
            void Go(const uint32_t ref);

            // Hot-path methods inlined here so the compiler sees them in every translation unit and can eliminate the call/return overhead.

            // Advance one character. The fast path (still inside the current in-memory buffer) is a single increment; load() is called only when
            // the buffer boundary is crossed - that happens rarely.
            void Forward() {
                if (++_pointer == _currentSize && _currentSize != 0) {
                    _startRef += _currentSize;
                    load(); // slow path: fetch next buffer chunk
                }
            }

            // Step back one character.  Fast path is a simple decrement.
            void Back() {
                if (_pointer > 0) { --_pointer; }
                else if (_startRef > 0) { --_startRef; load(); }
            }

            size_t GetRef()  const { return _startRef + _pointer; }
            char   Next();
            char   Prev()    const { return (_pointer > 0) ? _buffer[_pointer - 1] : prev; }
            bool   IsAtEnd() const { return _pointer == _currentSize && _currentSize == 0; }
            char   Get()     const { return _buffer[_pointer]; }

            void   Clear();
            size_t GetNoLines();
            IfcFileStream * Clone();
          private:
            void load();
            std::function<uint32_t(char *, size_t, size_t)> _dataSource;
            size_t _pointer     = 0;
            size_t _size;
            char   prev         = 0;
            size_t _currentSize = 0;
            size_t _startRef    = 0;
            char * _buffer;
            size_t noLines      = 0;
        };
        class IfcTokenChunk
        {
            public:
            	IfcTokenChunk(const size_t chunkSize, const size_t startRef, const size_t fileStartRef, IfcTokenStream::IfcFileStream *_fileStream);
              bool Clear(bool force);
              bool Clear();
              bool IsLoaded();
              size_t TokenSize();
              size_t GetTokenRef();
              void Push(void *v, const size_t size);
              size_t GetMaxSize();
              std::string_view ReadString(const size_t ptr,const size_t size); 
              template <typename T> T Read(const size_t ptr)
              {
                if (!_loaded) Load();
                T v;
                std::memcpy(&v, _chunkData+ptr, sizeof(T));
                return v;
              }
              template <typename T> void Push(T input)
              {
                Push(&input,sizeof(T));
              }
            private:
              void Load();
              bool _loaded=false;
              size_t _currentSize=0;
              size_t _startRef=0;
              size_t _fileStartRef;
              size_t _chunkSize;
            	uint8_t *_chunkData;
              IfcFileStream *_fileStream;
        };
        IfcTokenStream(size_t activeChunks, uint64_t maxChunks, std::vector<IfcTokenChunk> &chunks,IfcFileStream * fileStream);
        std::vector<IfcTokenChunk> _chunks;
        IfcTokenChunk * _cChunk;
        IfcFileStream * _fileStream;
  };
  
}