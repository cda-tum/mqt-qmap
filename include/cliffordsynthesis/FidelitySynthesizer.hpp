#include "cliffordsynthesis/CliffordSynthesizer.hpp"

#ifndef CS_FIDELITYSYNTHESIZER_HPP
    #define CS_FIDELITYSYNTHESIZER_HPP

class FidelitySynthesizer: public CliffordSynthesizer {
protected:
    void makeSpecificEncoding(const CliffordSynthesizer::SynthesisData& data, const SynthesisConfiguration& configuration) override;
    void updateResults(SynthesisResults& results) override;
};

#endif //CS_FIDELITYSYNTHESIZER_HPP
