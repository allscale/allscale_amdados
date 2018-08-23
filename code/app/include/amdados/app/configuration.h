//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017-2018
//-----------------------------------------------------------------------------

#pragma once

namespace amdados {

//=============================================================================
// Simple class for managing the application parameters.
// Parameters can be read from a text file where the format of each parameter
// is as follows:
//          <parameter name>  <parameter value>  # comment
// Parameter name is a single word without spaces, special symbols or
// quotation. Parameter value is a double, an integer or a string without
// spaces and quotation. Everything to the right of the symbol '#' is ignored.
// Blank lines are ignored as well.
//=============================================================================
class Configuration
{
private:
    enum ParamType { UNKNOWN_T, NUMERIC_T, STRING_T };

    // Simple placeholder of a parameter of any type.
    struct Parameter
    {
        std::string name;       ///< name of this parameter
        std::string svalue;     ///< "string" value, if applicable
        double      dvalue;     ///< "double" value, if applicable
        ParamType   type;       ///< type of this parameter

        Parameter() : name(), svalue(), dvalue(0.0), type(UNKNOWN_T) {}
    };

    std::map<std::string,Parameter> m_params;   ///< parameter placeholder

    void CheckExist(bool do_exist, const char * param_name) const;
    void CheckType(const Parameter & p, ParamType expected) const;
    void CheckIntRange(double value) const;

public:
    Configuration();
    void ReadConfigFile(const std::string & filename);
    void PrintParameters() const;

    int                 asInt(const char * param_name) const;
    size_t              asUInt(const char * param_name) const;
    double              asDouble(const char * param_name) const;
    const std::string & asString(const char * param_name) const;
    const char *        asCString(const char * param_name) const;

    bool IsExist(const char * param_name) const;
    bool IsInteger(const char * param_name) const;

    void SetInt(const char * param_name, int value = 0);
    void SetDouble(const char * param_name, double value = 0.0);
    void SetString(const char * param_name, const char * value = "");

	friend std::ostream& operator<<(std::ostream & out,
                                    const Configuration & c) {
		out << "Configuration: [ ";
		for (const auto & e : c.m_params) {
			out << e.first;
			out << ", ";
			out << e.second.name << ", ";
			out << e.second.svalue << ", ";
			out << e.second.dvalue << ", ";
			out << e.second.type << "; ";
		}
		out << "]";
		out << std::endl;
		return out;
	}
};

} // namespace amdados

