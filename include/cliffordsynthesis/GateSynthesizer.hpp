#include "cliffordsynthesis/CliffordSynthesizer.hpp"

#ifndef CS_GATESYNTHESIZER_HPP
    #define CS_GATESYNTHESIZER_HPP
class GateSynthesizer: public CliffordSynthesizer {
protected:
    void makeSpecificEncoding(const CliffordSynthesizer::SynthesisData& data, const SynthesisConfiguration& configuration) override;
    void updateResults(SynthesisResults& results) override;
};
#endif //CS_GATESYNTHESIZER_HPP
