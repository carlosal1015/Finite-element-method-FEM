#ifndef FREQUENCYSOURCE_HH
#define FREQUENCYSOURCE_HH

#include <iostream>
#include <cctype>
#include <iterator>
#include <algorithm>
#include <utility>

/**
 * Source of letters from a stream.
 *
 * Use it like this:
 *
 * auto source = streamLetterSource(std::cin);
 *
 * while (true)
 * {
 *   auto data = source.next();
 *   if (not data.second)
 *     break;
 *   // work with letter in data.first
 * }
 */
template<typename Stream>
class StreamLetterSource
{

public:

  // The type of data returned from next
  // Item::second tells you whether the data is valid,
  // Item::first contains the actual data if it is valid
  using Item = std::pair<char,bool>;

  StreamLetterSource(Stream& stream)
    : _stream(stream)
  {}

  // Get the next data item from the source.
  Item next()
  {
    unsigned char c;
    _stream >> c;
    return Item(c,bool(_stream));
  }

private:

  Stream& _stream;

};


/**
 * Construct a StreamLetterSource for a givn stream
 */
template<typename Stream>
StreamLetterSource<Stream> streamLetterSource(Stream& stream)
{
  return {{stream}};
}


/**
 * Source of words from a stream.
 *
 * Use it like this:
 *
 * auto source = streamWordSource(std::cin);
 *
 * while (true)
 * {
 *   auto data = source.next();
 *   if (not data.second)
 *     break;
 *   // work with word in data.first
 * }
 */
template<typename Stream>
class StreamWordSource
{

public:

  // The type of data returned from next
  // Item::second tells you whether the data is valid,
  // Item::first contains the actual data if it is valid
  using Item = std::pair<std::string,bool>;

  StreamWordSource(Stream& stream)
    : _stream(stream)
  {}

  // Get the next data item from the source.
  Item next()
  {
    std::string s;
    _stream >> s;
    if (not bool(_stream))
      return Item("",false);

    std::string out;
    std::copy_if(begin(s),end(s),back_inserter(out),[](auto c)
                 {
                   return std::isalpha(c);
                 });
    std::transform(begin(out),end(out),begin(out),[](auto c)
                   {
                     return std::tolower(c);
                   });
    return Item(std::move(out),true);
  }

private:

  Stream& _stream;

};


/**
 * Construct a StreamLetterSource for a givn stream
 */
template<typename Stream>
StreamWordSource<Stream> streamWordSource(Stream& stream)
{
  return {{stream}};
}

#endif // FREQUENCYSOURCE_HH
