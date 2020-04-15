/* ========================================================================
 phonebook.cpp
 * example for illustrating how to manipulate structure.
 *
 * takes a (MxN) structure matrix which has first field as
 * character array(name), and second field as scalar double (phone number).
 * This function returns a new structure (1x1)containing following fields:
 * for character array input, it will be (MxN) cell array;
 * and for numeric double (noncomplex, scalar) input, it will be (MxN)
 * cell array where each field is numeric array of type double.
 *
 * Build : from MATLAB
 *         >> mex phonebook.cpp
 * Usage with example : from MATLAB
 *         >> friends(1).name = 'Jordan Robert';
 *         >> friends(1).phone = 3386;
 *         >> friends(2).name = 'Mary Smith';
 *         >> friends(2).phone = 3912;
 *         >> friends(3).name = 'Stacy Flora';
 *         >> friends(3).phone = 3238;
 *         >> friends(4).name = 'Harry Alpert';
 *         >> friends(4).phone = 3077;
 *         >> phonebook(friends)
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2017 The MathWorks, Inc.
 *======================================================================= */

#include "mex.hpp"
#include "mexAdapter.hpp"
#include <string>
#include <memory>

using namespace matlab::mex;
using namespace matlab::data;

class MexFunction : public Function {
private:
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr;
public:
    // Constructor for the class. 
    MexFunction()
    {
        matlabPtr = getEngine();
    }

    // Helper function to print output string on MATLAB command prompt. 
    void displayOnMATLAB(std::ostringstream stream);

    // Helper function to generate an error message from given string,
    // and display it over MATLAB command prompt.
     
    void displayError(std::string errorMessage);

    // Helper function to information about an empty field in the structure. 
    void emptyFieldInformation(std::string fieldName, size_t index);

    // Helper function to information about an invalid field in the structure. 
    void invalidFieldInformation(std::string fieldName, size_t index);


    // This is the gateway routine for the MEX-file. 
    void
        operator()(ArgumentList outputs, ArgumentList inputs);

    // Make sure that the passed structure has valid data. 
    void checkStructureElements(StructArray const& matlabStructArray);

    // This function makes sure that user has provided structure as input,
    // and is not expecting more than one output in results.
     
    void checkArguments(ArgumentList outputs, ArgumentList inputs);
};