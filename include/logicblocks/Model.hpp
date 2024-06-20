#pragma once

#include "Logic.hpp"
#include "LogicBlock.hpp"
#include "LogicTerm.hpp"

#include <cstdint>

namespace logicbase {
class Model {
protected:
  Result result = Result::NDEF;

public:
  Model() = default;
  explicit Model(Result res) : result(res) {}
  virtual ~Model() = default;

  virtual int getIntValue(const LogicTerm& a, LogicBlock* lb) = 0;
  virtual bool getBoolValue(const LogicTerm& a, LogicBlock* lb) = 0;
  virtual double getRealValue(const LogicTerm& a, LogicBlock* lb) = 0;
  virtual uint64_t getBitvectorValue(const LogicTerm& a, LogicBlock* lb) = 0;
};
} // namespace logicbase
