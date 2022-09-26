#include "cliffordsynthesis/CliffordSynthesizer.hpp"

#ifndef QMAP_FIDELITYSYNTHESIZER_HPP
#define QMAP_FIDELITYSYNTHESIZER_HPP

class FidelitySynthesizer: public CliffordSynthesizer {
protected:
    void makeSpecificEncoding(const CliffordSynthesizer::SynthesisData& data) override;
    void updateResults(SynthesisResults& results) override;
};

#endif //QMAP_FIDELITYSYNTHESIZER_HPP
