#include "cliffordsynthesis/CliffordSynthesizer.hpp"

#ifndef CS_DEPTHSYNTHESIZER_HPP
    #define CS_DEPTHSYNTHESIZER_HPP

class DepthSynthesizer: public CliffordSynthesizer {
protected:
    void makeSpecificEncoding(const CliffordSynthesizer::SynthesisData& data, const SynthesisConfiguration& configuration) override;
    void updateResults(SynthesisResults& results) override;
};
#endif //CS_DEPTHSYNTHESIZER_HPP
