#include "fileparse.hpp"

namespace GFlowSimulation {

  FileParse::FileParse() : currentHead(nullptr), level(0) {};

  HeadNode* FileParse::parseFile(string configFile) {
    // Create a parent head
    HeadNode *root = new HeadNode;
    currentHead = root;
    // Create filestream
    std::ifstream fin(configFile);
    if (fin.fail()) {
      cout << "File parse creator failed to open file [" << configFile << "].\nExiting.\n";
      return nullptr;
    }
    // Do the actual parsing
    message = "Starting parse...\n";
    level = 0;
    getBody(fin);
    // Set current head back to null
    currentHead = nullptr;
    // Return the head
    return root;
  }

  inline void FileParse::passComment(std::ifstream& fin, bool mline) {
    string comment;
    // Start right after "//" or "/*"
    char c;
    fin.get(c);
    while (!fin.eof()) {
      if (mline && c=='*') {
        fin.get(c);
        if (fin.eof()) {
          // Message
          message += (tabs() + (mline ? "ML" : "SL") + " <" + comment +">\n");
          message += "--> EOF\n";
          // End of file encountered
          return; // Check end of file
        }
        if (c=='/') {
          // Message
          message += (tabs() + (mline ? "ML" : "SL") + " <" + comment +">\n");
          // End of the comment
          return; 
        }
      }
      else if (c=='\n' || c=='\r') {// Single line comments end with at newline
        if (mline) comment += "[\\n]";
        else {
          // Message
          message += (tabs() + (mline ? "ML" : "SL") + " <" + comment +">\n");
          // In case a newline signifies something for whoever called this function
          fin.putback(c); 
          return;
        }
      }
      else comment.push_back(c);
      // Get next character
      fin.get(c);
    }
  }

  inline bool FileParse::passSpaces(std::ifstream& fin) {
    char c;
    // Get first character
    fin.get(c);
    // Loop
    while (!fin.eof()) {
      if (c==' ');
      else if (c=='\n' || c=='\r') return true;
      else { // Encountered a non-whitespace
        fin.putback(c);
        return false;
      }
      // Get next character
      fin.get(c);
    }
  }

  inline void FileParse::passWhiteSpaces(std::ifstream& fin) {
    char c;
    // Get first character
    fin.get(c);
    // Loop
    while (!fin.eof()) {
      if (c==' ' || c=='\n' || c=='\r');
      else { // Encountered a non-whitespace
        fin.putback(c);
        return;
      }
      // Get next character
      fin.get(c);
    }
  }

  inline void FileParse::getBody(std::ifstream& fin) {
    // Look for heads
    char c;
    bool end = false;
    while (!fin.eof() && !end) {
      passWhiteSpaces(fin);

      if (!fin.eof()) fin.get(c);
      else return;

      if (c=='}') // End of a body
        return;
      else if (c=='/') { // Could be the start of a comment
        checkComment(fin);
      }
      else {
        fin.putback(c);
        getHead(fin);
      }
    }
  }

  inline void FileParse::getHead(std::ifstream& fin) {
    if (fin.eof()) return;
    // Create a new head node
    HeadNode *node = new HeadNode;

    // Look for the heading. Ends with ":", or "{". If "{" is encountered, body = true.
    bool body = false;
    getHeading(fin, node->heading, body);

    // Message
    message += (tabs() + "Level " + toStr(level) + ": Heading: [" + node->heading + "]\n");

    // Set node as the current head
    node->parent = currentHead;
    currentHead = node;

    // Look for parameters
    bool newLine = passSpaces(fin);
    
    // Get the parameters - adds them to the current head node
    if (!fin.eof() && !newLine && !body)
      body = getParameters(fin);

    ++level;
    if (body && !fin.eof()) getBody(fin);
    --level;

    // Add this node to the parent node
    node->parent->subHeads.push_back(node);

    // Return to parent node
    currentHead = node->parent;
  }

  inline void FileParse::getHeading(std::ifstream& fin, string& heading, bool& found) {
    char c;
    fin.get(c);
    bool whitespace = false;
    while (!fin.eof()) {
      if (c=='{') {
        found = true;
        return;
      }
      else if (c==':') {
        found = false;
        return;
      }
      else if (c=='\n' || c=='\r' || c==' ') whitespace = true;
      else { // Just a regular character. Could be a whitespace
        if (whitespace) { // Error
          throw UnexpectedToken();
        }
        else heading.push_back(c);
      }
      // Get the next character
      fin.get(c);
    }
  }

  inline bool FileParse::getParameters(std::ifstream& fin) {    
    bool more = true, body = false;
    ++level;
    while (more) getParam(fin, more, body);
    --level;
    // Return whether we expect a body or not
    return body;
  }

  inline void FileParse::getParam(std::ifstream& fin, bool& more, bool& body) {
    char c;
    bool end = false, a_part = true;
    fin.get(c);
    string a(""), b("");
    while (!fin.eof()) {
      if (c=='/') {
        checkComment(fin);
      }
      else if (c=='=') {
        a_part = false;
      }
      else if (c==',') {
        more = true;
        break;
      }
      else if (c=='{') {
        more = false;
        body = true;
        break;
      }
      else if (c=='\n' || c=='\r') {
        body = false;
        more = false;
        break;
      }
      else if (c==' ');
      else if (c=='\n' || c=='\r') {
        more = false;
        body = false;
        break;
      }
      else {
        if (a_part) a.push_back(c);
        else        b.push_back(c);
      }

      // Get next character
      if (!fin.eof()) fin.get(c);
    }
    // Set param
    if (!a.empty() || !b.empty()) {
      // Message
      message += (tabs() + "*> Param: [" + a + "]");
      if (!b.empty()) message += (", [" + b + "]");
      message += "\n";
      // Push nodes
      currentHead->params.push_back(new ParamNode(a, b));
    }
  }

  inline void FileParse::checkComment(std::ifstream& fin) {
    char c;
    fin.get(c);
    // Make sure we are not at eof
    if (fin.eof()) return;
    // Check what kind of comment this is
    if (c=='/')      passComment(fin, false);
    else if (c=='*') passComment(fin, true);
    else             throw UnexpectedToken();
  }

  inline string FileParse::tabs() {
    string tbs = "";
    for (int i=0; i<level; ++i) tbs += "\t";
    return tbs;
  }

}