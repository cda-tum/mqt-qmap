#include "cliffordsynthesis/CliffordSynthesizer.hpp"

#ifndef QMAP_FIDELITYSYNTHESIZER_HPP
#define QMAP_FIDELITYSYNTHESIZER_HPP

class Fidelitysynthesizer : public CliffordSynthesizer {
protected:
    void makeSynthesis(const CliffordSynthesizer::SynthesisData& data) override;
};

#endif //QMAP_FIDELITYSYNTHESIZER_HPP
