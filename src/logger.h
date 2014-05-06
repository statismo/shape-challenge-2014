#ifndef __LOGGER_H
#define __LOGGER_H

#include <fstream>

// Extremely simplistic logging class
enum LogLevel {logDEBUG, logINFO, logERROR};

class Logger { 
public:
	Logger(const char* filename, LogLevel logLevel) : m_logLevel(logLevel) { 
		m_os.open(filename);
		if (!m_os) { 
			std::ostringstream omsg;
			omsg << "could not open logging file with name " << filename;
			throw std::exception(omsg.str().c_str());
		}
	}

	~Logger() { m_os.close(); } 

	std::ofstream& Get(const LogLevel& logLevel) { 
		// we currently ignore the logLevel
		return m_os; 
	}  

private:
	std::ofstream m_os;
	LogLevel m_logLevel;
};

#endif