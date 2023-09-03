#pragma once

#ifndef ERROR_CODES_H
#define ERROR_CODES_H

#include <string>

/* Error codes */
typedef enum
{
    ERR_NoError = 0,                // No error
    ERR_BadInputVecs = -6201,       // Corrupted input vectors
    ERR_ZeroInputVals = -6202,      // Input data error, some values equals 0
    ERR_BadCutOff = -6203,          // Ñut-off points must be less than 0.9 in total and more than 0 each
    ERR_BadLinFit = -6204,          // The length of linear fit must be less than 0.9
    ERR_BadFactorizing = -6205,     // Error after Ramp|Signal factorizing
    ERR_BadNoise = -6206,           // Error after noise extracting
    ERR_BadSegInput = -6207,        // Input segment's values error while segmend approximating
    ERR_TooFewSegs = -6208,         // Less then 4 segments found, check input arrays  
    ERR_BadSegsLength = -6209,      // Error in finding segments length, check input params
    ERR_BadLinearRamp = -6210,      // Error in Ramp linearizing, check cut-off params
    ERR_TooManyAttempts = -6211,    // More than 5 attempts to find signal, check if signal is noise or not
    ERR_BadStartEnd = -6212,        // Error in finding start|end of signal, check if signal is noise or not
    ERR_TooManySegs = -6213,        // Too many segments, check if signal is noise or not
    ERR_NoSegs = -6214,             // No segments found, check if signal is noise or not
    ERR_BadDiagNum = -6215,         // Diagnostics number must be > 0 and < 2
    ERR_IdxOutOfRange = -6216,      // Index is out of range. Programm continued running
    ERR_Exception = -6217,          // An exception been occured while running, script did not stopt working
    ERR_BadStEndTime = -6218,       // Error!: start time must be less then end time and total time, more than 0\n\
                                            end time must be less then total time, more then 0
    ERR_BufferExtension = -6219,    // Buffer extension, bad impulse
    ERR_StorageIsEmpty = -6220,     // No data been set to work with. Check input arrays
    ERR_BadExpression = -6221,      // Error at passed expression
    ERR_ExpressionNotSet = -6222,   // Setting up of initial values of variables is only allowed after setting up the expression
    ERR_NoetAllDataSet = -6223      // Running calculation is only allowed after setting up voltage data, signal data and additional data
} error_codes;

/* return the description of a particular error code */
inline std::string ERR_GetErrorDescription(int err)
{
    switch (err)
    {
    case ERR_BadInputVecs:
        return "Corrupted input vectors";
    case ERR_ZeroInputVals:
        return "Input data error, some values equals 0";
    case ERR_BadCutOff:
        return "Ñut-off points must be less than 0.9 in total and more than 0 each";
    case ERR_BadLinFit:
        return "The length of linear fit must be less than 0.9";
    case ERR_BadFactorizing:
        return "Error after Ramp|Signal factorizing";
    case ERR_BadNoise:
        return "Error after noise extracting";
    case ERR_BadSegInput:
        return "Input segment's values error while segmend approximating";
    case ERR_TooFewSegs:
        return "Less then 4 segments found, check input arrays";
    case ERR_BadSegsLength:
        return "Error in finding segments length, check input params";
    case ERR_BadLinearRamp:
        return "Error in Ramp linearizing, check cut-off params";
    case ERR_TooManyAttempts:
        return "More than 5 attempts to find signal, check if signal is noise or not";
    case ERR_BadStartEnd:
        return "Error in finding start|end of signal, check if signal is noise or not";
    case ERR_TooManySegs:
        return "Too many segments, check if signal is noise or not";
    case ERR_NoSegs:
        return "No segments found, check if signal is noise or not";
    case ERR_BadDiagNum:
        return "Diagnostics number must be > 0 and < 3";
    case ERR_IdxOutOfRange:
        return "Index is out of range. Programm continued running";
    case ERR_Exception:
        return "An exception been occured while running, script did not stoped working";
    case ERR_BadStEndTime:
        return "Error!: start time must be less then end time and total time, more than 0\n\
end time must be less then total time, more then 0";
    case ERR_BufferExtension:
        return "Buffer extension, bad impulse";
    case ERR_StorageIsEmpty:
        return "No data been set to work with. Check input arrays";
    case ERR_ExpressionNotSet:
        return "Setting up of initial values of variables is only allowed after setting up the expression";
    case ERR_NoetAllDataSet:
        return "Running calculation is only allowed after setting up voltage data, signal data and additional data";
    default:
        return "No Error";
    }
};

/* std::string input overload */
inline std::string ERR_GetErrorDescription(std::string s_err)
{
    return s_err;
};

inline int err_num(int err)
{
    return err;
};

inline int err_num(std::string s_err)
{
    return -1;
};

#endif
