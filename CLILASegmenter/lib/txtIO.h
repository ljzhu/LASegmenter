#ifndef txtIO_h_
#define txtIO_h_

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

namespace douher
{
  /*============================================================
   * readTextLineToListOfString   
   */
  std::vector< std::string > readTextLineToListOfString(const char* textFileName)
  {
    /* The file MUST end with an empty line, then each line will be
       stored as an element in the returned vector object. */


    // Here speed and space is not a big issue, so just let it copy and return stuff.
    std::vector< std::string > listOfStrings;

    std::ifstream f(textFileName);

    int maxFileNameLength = 1024;
    char n[maxFileNameLength];

    if (f.good())
      {
        while(!f.eof())
          {
            f.getline(n, maxFileNameLength);
            std::string name(n);
            listOfStrings.push_back(name);
          }
      }


    // the last line is an empty line, so remove
    listOfStrings.pop_back();

    f.close();

    return listOfStrings;
  }



  /*============================================================
   * readNumbersToVector
   */
  void readNumbersToVector(const char* textFileName, std::vector< double >& v)
  {
    /* The file MUST end with an empty line, then each line will be
       stored as an element in the returned vector object. */


    // v will be cleared !!!
    v.clear();

    std::ifstream f(textFileName);

    if (f.good())
      {
        while(!f.eof())
          {
            double a;
            f >> a;
            v.push_back(a);
          }
      }


    // the last line is an empty line, so remove
    v.pop_back();

    f.close();

    return;
  }



} //douher


#endif
