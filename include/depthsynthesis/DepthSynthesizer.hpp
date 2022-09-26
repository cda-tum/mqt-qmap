#include "cliffordsynthesis/CliffordSynthesizer.hpp"

#ifndef QMAP_DEPTHSYNTHESIZER_HPP
    #define QMAP_DEPTHSYNTHESIZER_HPP

class DepthSynthesizer: public CliffordSynthesizer {
protected:
    void makeSpecificEncoding(const CliffordSynthesizer::SynthesisData& data, const SynthesisConfiguration& configuration) override;
    void updateResults(SynthesisResults& results) override;
};
#endif //QMAP_DEPTHSYNTHESIZER_HPP
