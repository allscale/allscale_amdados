//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017-2018
//-----------------------------------------------------------------------------

#include "allscale/utils/assert.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <limits>
#include <map>
#include "../include/debugging.h"
#include "../include/configuration.h"

namespace amdados {

//-----------------------------------------------------------------------------
// Constructor.
//-----------------------------------------------------------------------------
Configuration::Configuration() : m_params()
{
}

//-----------------------------------------------------------------------------
// Function reads parameters from the configuration file.
//-----------------------------------------------------------------------------
void Configuration::ReadConfigFile(const std::string & filename)
{
    std::fstream f(filename, std::ios::in);
    assert_true(f.good()) << "ERROR: failed to open configuration file"
                          << std::endl;
    int lineNo = 1;
    std::string line, token, name;
    // Read line by line.
    while (std::getline(f, line)) {
        if (line.empty())
            continue;
        std::stringstream ss(line);
        // Parse the line.
        for (int count = 1; ss >> token; ++count) {
            assert_true(count <= 3)
                << "ERROR at line " << filename << ":" << lineNo << std::endl
                << "the valid line layout: "
                << "'<comment>' or '<name> <value> <comment>'" << std::endl;

            // Comment is the 1st or the 3rd token in a line.
            if (!token.empty() && (token[0] == '#')) {
                assert_true((count == 1) || (count == 3))
                    << "ERROR at line " << filename << ":" << lineNo
                    << std::endl << "comment may present in the "
                    << "1st or the 3rd token in a line" << std::endl;
                break;
            }
            // If not comment, the 1st token is the name, the 2nd one is
            // the value. The value is a "double" if sscanf() succeeded,
            // otherwise a "string".
            if (count == 1) {
                name = token;
            } else if (count == 2) {
                double value = 0.0;
                if (std::sscanf(token.c_str(), "%lf", &value) == 1) {
                    SetDouble(name.c_str(), value);
                } else {
                    SetString(name.c_str(), token.c_str());
                }
            }
        }
        assert_true(!f.bad())
            << "ERROR at line " << filename << ":" << lineNo << std::endl
            << "failure while reading configuration file" << std::endl;
    }
}

//-----------------------------------------------------------------------------
// Function prints out all the parameters in debugging mode,
// otherwise prints nothing.
//-----------------------------------------------------------------------------
void Configuration::PrintParameters() const
{
#ifdef AMDADOS_DEBUGGING
    std::cout << std::endl << "Parameters:" << std::endl;
    for (const auto & p : m_params) {
        std::cout << p.second.name << " : ";
        switch (p.second.type) {
            case NUMERIC_T : std::cout << p.second.dvalue << std::endl;
                             break;
            case STRING_T  : std::cout << p.second.svalue << std::endl;
                             break;
            default        : assert_fail() << "never get here" << std::endl;
                             break;
        }
    }
    std::cout << std::endl;
#endif
}

//-----------------------------------------------------------------------------
// Function checks parameter existence flag.
//-----------------------------------------------------------------------------
void Configuration::CheckExist(bool do_exist, const char * param_name) const
{
    assert_true(do_exist) << "Parameter " << param_name
                          << " was not found" << std::endl;
}

//-----------------------------------------------------------------------------
// Function checks parameter's type.
//-----------------------------------------------------------------------------
void Configuration::CheckType(const Configuration::Parameter & p,
                              ParamType expected) const
{
    assert_true(p.type != UNKNOWN_T)
        << "Parameter " << p.name << " has no type defined" << std::endl;
    assert_true(p.type == expected)
        << "Parameter " << p.name << " is expected to be a "
        << ((expected == NUMERIC_T) ? "numeric" : "string") << std::endl;
}

//-----------------------------------------------------------------------------
// Function checks parameter's value does fit the integer range.
//-----------------------------------------------------------------------------
void Configuration::CheckIntRange(double value) const
{
    assert_true((std::numeric_limits<int>::lowest() <= value) &&
                (value <= std::numeric_limits<int>::max()))
        << "the value '" << value << "' causes numerical overflow" << std::endl
        << "when converted to an integer parameter" << std::endl;
}

//-----------------------------------------------------------------------------
// Function returns the integer parameter value.
//-----------------------------------------------------------------------------
int Configuration::asInt(const char * param_name) const
{
    auto it = m_params.find(param_name);
    CheckExist(it != m_params.end(), param_name);
    CheckType(it->second, NUMERIC_T);
    CheckIntRange(it->second.dvalue);
    return static_cast<int>(std::floor(it->second.dvalue + 0.5));
}

//-----------------------------------------------------------------------------
// Function returns the unsigned integer parameter value.
//-----------------------------------------------------------------------------
size_t Configuration::asUInt(const char * param_name) const
{
    int v = asInt(param_name);
    assert_true(v >= 0);
    static_assert(sizeof(int) <= sizeof(size_t), "wrong size of built-ins");
    return static_cast<size_t>(v);
}

//-----------------------------------------------------------------------------
// Function returns the floating-point parameter value.
//-----------------------------------------------------------------------------
double Configuration::asDouble(const char * param_name) const
{
    auto it = m_params.find(param_name);
    CheckExist(it != m_params.end(), param_name);
    CheckType(it->second, NUMERIC_T);
    return it->second.dvalue;
}

//-----------------------------------------------------------------------------
// Function returns the string parameter value.
//-----------------------------------------------------------------------------
const std::string & Configuration::asString(const char * param_name) const
{
    auto it = m_params.find(param_name);
    CheckExist(it != m_params.end(), param_name);
    CheckType(it->second, STRING_T);
    return it->second.svalue;
}

//-----------------------------------------------------------------------------
// Function returns the string-pointer parameter value.
//-----------------------------------------------------------------------------
const char * Configuration::asCString(const char * param_name) const
{
    auto it = m_params.find(param_name);
    CheckExist(it != m_params.end(), param_name);
    CheckType(it->second, STRING_T);
    return it->second.svalue.c_str();
}

//-----------------------------------------------------------------------------
// Function returns "true" if parameter has an integer value.
//-----------------------------------------------------------------------------
bool Configuration::IsInteger(const char * param_name) const
{
    return (asDouble(param_name) == asInt(param_name));
}

//-----------------------------------------------------------------------------
// Function initializes an integer parameter.
//-----------------------------------------------------------------------------
void Configuration::SetInt(const char * param_name, int value)
{
    SetDouble(param_name, static_cast<double>(value));
    CheckIntRange(value);
}

//-----------------------------------------------------------------------------
// Function initializes a floating-point parameter.
//-----------------------------------------------------------------------------
void Configuration::SetDouble(const char * param_name, double value)
{
    auto it = m_params.find(param_name);
    if (it == m_params.end()) {
        Parameter & p = m_params[param_name];
        p.name = param_name;
        p.dvalue = value;
        p.type = NUMERIC_T;
    } else {
        CheckType(it->second, NUMERIC_T);
        it->second.dvalue = value;
    }
}

//-----------------------------------------------------------------------------
// Function initializes a string parameter.
//-----------------------------------------------------------------------------
void Configuration::SetString(const char * param_name, const char * value)
{
    auto it = m_params.find(param_name);
    if (it == m_params.end()) {
        Parameter & p = m_params[param_name];
        p.name = param_name;
        p.svalue = value;
        p.type = STRING_T;
    } else {
        CheckType(it->second, STRING_T);
        it->second.svalue = value;
    }
}

} // namespace amdados

