
#ifndef INCLUDE_XML_HELPER_H_
#define INCLUDE_XML_HELPER_H_

#include <constants.h>

class XPathException: public std::exception
{
  virtual const char* what() const throw()
  {
    return "No match for XPath expression.";
  }
} xPathException;

int loadXml(std::string& fxml, pugi::xml_document& doc)
{
	pugi::xml_parse_result result = doc.load_file(fxml.c_str());
	if (!result)
	{
		std::cout << "\nError:  cannot open or parse GPX file \""
				<< fxml << "\"." << std::endl;
		return ERR_XML_OPEN;
	}
	return 0;
}

#endif  // INCLUDE_XML_HELPER_H_

