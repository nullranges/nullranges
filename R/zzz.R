# suppress R CMD CHECK note on globals
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))
