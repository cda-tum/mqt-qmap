#include "cliffordsynthesis/CliffordSynthesizer.hpp"

#ifndef QMAP_GATESYNTHESIZER_HPP
#define QMAP_GATESYNTHESIZER_HPP
class GateSynthesizer: public CliffordSynthesizer {
protected:
    void makeSynthesis(const CliffordSynthesizer::SynthesisData& data) override;
};
#endif //QMAP_GATESYNTHESIZER_HPP
