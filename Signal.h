#ifndef Signal_h
#define Signal_h

#include "Rtypes.h"
#include <iostream>

struct Signal
{
  int Label;
  double Time;
  int Channel;
  bool isValid;

  Signal() : isValid(false) {}
  Signal(int label, double time, int channel) : Label(label), Time(time), Channel(channel), isValid(true) {}
  bool operator==(const Signal &other) const
  {
    return Label == other.Label && Time == other.Time && Channel == other.Channel;
  }
  bool operator!=(const Signal &other) const
  {
    return Label != other.Label || Time != other.Time || Channel != other.Channel;
  }

  friend std::ostream &operator<<(std::ostream &os, const Signal &s)
  {
    os << "Label:" << s.Label << "\t Time: " << s.Time << "\t Channel: " << s.Channel;
    return os;
  }

  ClassDef(Signal, 1)  // Signal
};

#endif