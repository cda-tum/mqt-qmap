#include "cliffordsynthesis/CliffordSynthesizer.hpp"

#ifndef QMAP_DEPTHSYNTHESIZER_HPP
    #define QMAP_DEPTHSYNTHESIZER_HPP

class DepthSynthesizer: public CliffordSynthesizer {
protected:
    void makeSynthesis(const CliffordSynthesizer::SynthesisData& data) override;
};
#endif //QMAP_DEPTHSYNTHESIZER_HPP
