#ifndef Signal_h
#define Signal_h

#include "Rtypes.h"
#include <iostream>

struct Signal
{
  int Label;
  double Time;
  double Channel;
  int Multiplicity = 1;
  bool isValid;

  Signal() : Label(-1), Time(-1), Channel(0), Multiplicity(0), isValid(false) {}
  Signal(int label, double time, int channel) : Label(label), Time(time), Channel(channel), Multiplicity(1), isValid(true) {}
  Signal(int label, double time, double channel, int multiplicity) : Label(label), Time(time), Channel(channel), Multiplicity(multiplicity), isValid(true) {}
  bool operator==(const Signal &other) const
  {
    return Label == other.Label && Time == other.Time && Channel == other.Channel;
  }
  bool operator!=(const Signal &other) const
  {
    return Label != other.Label || Time != other.Time || Channel != other.Channel;
  }

  Signal operator+(const Signal &other)
  {
    return Signal(Label, Time, Channel + other.Channel, Multiplicity + other.Multiplicity);
  }

  Signal operator+=(const Signal &other)
  {
    if (isValid && other.isValid)
    {
      Channel += other.Channel;
      Multiplicity += other.Multiplicity;
      return *this;
    }
    else if (isValid && !other.isValid)
    {
      return *this;
    }
    else if (!isValid && other.isValid)
    {
      Label = other.Label;
      Time = other.Time;
      Channel = other.Channel;
      Multiplicity = other.Multiplicity;
      isValid = true;
      return *this;
    }
    else
    {
      return Signal();
    }
  }

  Signal operator/(const double &value)
  {
    return Signal(Label, Time, Channel / value, Multiplicity);
  }

  friend Signal operator/(double value, Signal &signal)
  {
    return signal / value;
  }

  Signal operator*(const double &value)
  {
    return Signal(Label, Time, Channel * value, Multiplicity);
  }

  friend Signal operator*(double value, Signal &signal)
  {
    return signal * value;
  }

  friend std::ostream &operator<<(std::ostream &os, const Signal &s)
  {
    os << "Label:" << s.Label << "\t Time: " << s.Time << "\t Channel: " << s.Channel;
    return os;
  }

  ClassDef(Signal, 1) // Signal
};

#endif