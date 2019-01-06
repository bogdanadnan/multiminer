#include "testcase.h"

static char toHex(std::uint8_t digit) {
    return digit >= 10 ? 'a' + (digit - 10) : '0' + digit;
}

static void dumpBytes(std::ostream &out, const void *data, std::size_t size)
{
    auto bdata = static_cast<const std::uint8_t *>(data);
    for (std::size_t i = 0; i < size; i++) {
        auto val = bdata[i];
        out << toHex((val >> 4) & 0xf) << toHex(val & 0xf);
    }
}

void TestCase::dump(std::ostream &out) const
{
    out << "t=" << params.getTimeCost()
        << " m=" << params.getMemoryCost()
        << " p=" << params.getLanes()
        << " pass=";
    dumpBytes(out, input, inputLength);

    if (params.getSaltLength()) {
        out << " salt=";
        dumpBytes(out, params.getSalt(), params.getSaltLength());
    }

    if (params.getAssocDataLength()) {
        out << " ad=";
        dumpBytes(out, params.getAssocData(), params.getAssocDataLength());
    }

    if (params.getSecretLength()) {
        out << " secret=";
        dumpBytes(out, params.getSecret(), params.getSecretLength());
    }
}
