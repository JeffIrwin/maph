
#ifndef INCLUDE_XML_HELPER_H_
#define INCLUDE_XML_HELPER_H_

class XPathException: public std::exception
{
  virtual const char* what() const throw()
  {
    return "No match for XPath expression.";
  }
} xPathException;

#endif  // INCLUDE_XML_HELPER_H_

