#include <iostream>
#include <fstream>
#include <string>

using namespace std;

struct EventId {
    /* Stores the information needed to distinguish events.
     * The < operator is usd to sort events, the actual order is not of great importance.
     */
    EventId( unsigned int eventNumber_, unsigned int lumiBlockNumber_, unsigned int runNumber_ ) :
        eventNumber(eventNumber_),
        lumiBlockNumber(lumiBlockNumber_),
        runNumber(runNumber_) {}
    unsigned int eventNumber, lumiBlockNumber, runNumber;
    bool operator < ( const EventId& rh ) const {
        return ( eventNumber != rh.eventNumber ? eventNumber < rh.eventNumber :
                (lumiBlockNumber != rh.lumiBlockNumber ? lumiBlockNumber < rh.lumiBlockNumber :
                 runNumber < rh.runNumber ) );
    }
};

class EventIdCheck {
  public:
    EventIdCheck() {}
    EventIdCheck(const string& filename) {
      ifstream myfile(filename);
      string line;
      if (myfile.is_open()) {
        while (getline(myfile, line)) {
          smatch sm;
          regex e("(\\d+):(\\d+):(\\d+).*");
          regex_match(line, sm, e);
          if (sm.size() != 4) {
            cout << "Can not match line " << line << endl;
            continue;
          }
          unsigned run = stoul(sm[1]);
          unsigned lum = stoul(sm[2]);
          unsigned evt = stoul(sm[3]);
          if( !ids.insert( EventId( evt, lum, run ) ).second ) {
            cout << "Warning, event " << evt << " already in set" << endl;
          }
        }
      } else {
        cerr << "Unable to open file " << filename << endl;
      }
    }
    ~EventIdCheck() {}

    bool check(unsigned run, unsigned lum, unsigned evt) {
      return ids.find(EventId(evt, lum, run)) != ids.end();
    }
  private:
    set<EventId> ids;
};
