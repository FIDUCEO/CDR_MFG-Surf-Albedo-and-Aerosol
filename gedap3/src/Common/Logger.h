/**********************************************************************
* The copyrights for the GEDAP algorithm and computer codes, remain with 
* Rayference SPRL as an Intellectual Property Right 
*
* Any use, in source and binary forms, distribution, modification, for 
* commercial, business, research and any other purposes is
*
*                           *** STRICTLY FORBIDDEN ***
*
* GEDAP may not be distributed or sold to any other commercial, 
* business, research or other partners under any circumstances.
* In any case this comment is part of the code as Software Legal Information,
* it cannot be modified and must be kept as header for all source codes in 
* which it appears.
*
* All documents and software are protected under the international Copyright 
* Laws, and Rayference SPRL reserves all rights.
*
*
* Author       :     Rayference Copyright (c) 
*
*********************************************************************/
#ifndef LOGGER_H_
#define LOGGER_H_
#include <iostream>
#include <map>
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
namespace py = pybind11;
enum LogLevel {
    _CRITICAL_ = 50,
    _ERROR_ = 40,
    _WARNING_ = 30,
    _INFO_ = 20,
    _DEBUG_ = 10,
    _NOTSET_ = 0
};
static std::map<LogLevel, const char *> LogLevelString = {
    {_CRITICAL_, "CRITICAL"},
    {_ERROR_, "ERROR"},
    {_WARNING_, "WARNING"},
    {_INFO_, "INFO"},
    {_DEBUG_, "DEBUG"},
    {_NOTSET_, "NOTSET"}
};
class AbstractLogger {
public:
    AbstractLogger() = default;
    virtual ~AbstractLogger() = default;
    virtual void write_log(LogLevel level, std::string class_name, std::string method_name, std::string message) = 0;
};
class PythonLogger : public AbstractLogger {
private:
    py::object logger;
public:
    explicit PythonLogger(py::object logger) {
        this->logger = logger;
    }
    void write_log(LogLevel level, std::string class_name, std::string method_name, std::string message) override {
        logger(static_cast<int>(level), class_name, method_name, message);
    }
};
class CoutLogger : public AbstractLogger {
public:
    CoutLogger() = default;
    void write_log(LogLevel level, std::string class_name, std::string method_name, std::string message) override {
        std::cout << LogLevelString[level] << ":" << "cpp"
                  << " [" << class_name << "::" << method_name << "] "
                  << message
                  << std::endl;
    }
};
#endif
