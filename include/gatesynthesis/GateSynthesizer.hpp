#include "cliffordsynthesis/CliffordSynthesizer.hpp"

#ifndef QMAP_GATESYNTHESIZER_HPP
#define QMAP_GATESYNTHESIZER_HPP
class GateSynthesizer: public CliffordSynthesizer {
protected:
    void makeSpecificEncoding(const CliffordSynthesizer::SynthesisData& data) override;
    void updateResults(SynthesisResults& results) override;
};
#endif //QMAP_GATESYNTHESIZER_HPP
