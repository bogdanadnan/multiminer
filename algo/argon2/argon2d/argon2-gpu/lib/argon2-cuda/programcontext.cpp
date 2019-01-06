#include "programcontext.h"
#include "kernels.h"

#define THREADS_PER_LANE 32

namespace argon2 {
namespace cuda {

ProgramContext::ProgramContext(
        const GlobalContext *globalContext,
        const std::vector<Device> &devices,
        Type type, Version version)
    : globalContext(globalContext), type(type), version(version)
{
}

} // namespace cuda
} // namespace argon2

