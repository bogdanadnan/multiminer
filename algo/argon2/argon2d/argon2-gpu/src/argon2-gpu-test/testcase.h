#ifndef TESTCASE_H
#define TESTCASE_H

#include "argon2-gpu-common/argon2params.h"

#include <iostream>

class TestCase
{
private:
    argon2::Argon2Params params;
    const void *output;
    const void *input;
    std::size_t inputLength;

public:
    const argon2::Argon2Params &getParams() const { return params; }
    const void *getOutput() const { return output; }
    const void *getInput() const { return input; }
    std::size_t getInputLength() const { return inputLength; }

    TestCase(const argon2::Argon2Params &params, const void *output,
             const void *input, std::size_t inputLength)
        : params(params), output(output),
          input(input), inputLength(inputLength)
    {
    }

    void dump(std::ostream &out) const;
};

#endif // TESTCASE_H
