#ifndef __READLATTPAR__
#define __READLATTPAR__


// a slight modification on ReadLattice class (if not inheritance)
// to enable parallel reading/writing of "Gauge Connection Format" lattice data

#include <stdlib.h>	// exit()
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <map>
#include <string>

#include <util/gjp.h>
#include <util/vector.h>
#include <util/lattice.h>
#include <util/data_types.h>
#include <util/verbose.h>
#include <util/error.h>
#include <alg/do_arg.h>
#include <util/qioarg.h>
#include <util/fpconv.h>
#include <util/iostyle.h>
#include <util/latheader.h>

CPS_START_NAMESPACE
using namespace std;


class  ReadLatticeParallel : public QioControl{
  // which determines parallel reading or serial reading

 private:
  char *cname;
  LatticeHeader hd;

public:
  // ctor for 2-step loading
  ReadLatticeParallel()
    : QioControl(), cname("ReadLatticeParallel"), UseParIO(1)
    {  }

  // ctor invoking loading behavior
  ReadLatticeParallel(Lattice & lat, const char * filename, const Float chkprec = 0.01)
    : 
    QioControl(),
    cname("ReadLatticeParallel") , 
    UseParIO(1)
    {        
    QioArg rd_arg(filename,chkprec);
    read(lat,rd_arg);
  }

  // ctor invoking loading behavior
  ReadLatticeParallel(Lattice & lat, const QioArg & rd_arg) 
    : QioControl(), cname("ReadLatticeParallel"), UseParIO(1)
  {
    read(lat,rd_arg);
  }
  
  virtual ~ReadLatticeParallel() {}

  void read(Lattice & lat, const char * filename, const Float chkprec = 0.01){
    QioArg rd_arg(filename,chkprec);
    read(lat,rd_arg);
  }

  void read(Lattice & lat, const QioArg & rd_arg);

 private:
  bool CheckPlaqLinktrace(Lattice & lat, const QioArg & rd_arg,
			  const Float plaq_inheader, const Float linktrace_inheader);

 private:
    bool UseParIO;
 public:
    inline void setParallel() { UseParIO = 1; }
    inline void setSerial() { 
#if TARGET == QCDOC
      UseParIO = 0; 
#else
      const char * fname = "setSerial()";
      VRB.Flow(cname,fname,"On non-QCDOC platform, setSerial() has no effect!\n");
#endif
    }
    inline int parIO() const { return UseParIO; }
};


class ReadLatticeSerial : public ReadLatticeParallel {
 private:
  char * cname;

 public:
  // ctor for 2-step loading
  ReadLatticeSerial() : ReadLatticeParallel(), cname("ReadLatticeSerial") {
    setSerial();
  }

  // ctor invoking loading behavior
  ReadLatticeSerial(Lattice & lat, const char * filename, const Float chkprec = 0.01) 
    : ReadLatticeParallel(), cname("ReadLatticeSerial") {
    setSerial();
    read(lat, filename, chkprec);
  }

  // ctor invoking loading behavior
  ReadLatticeSerial(Lattice & lat, const QioArg & rd_arg) 
    : ReadLatticeParallel(), cname("ReadLatticeSerial") {
    setSerial();
    read(lat, rd_arg);
  }
  
  virtual ~ReadLatticeSerial() {}
};


CPS_END_NAMESPACE
#endif 


