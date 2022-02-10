use Time;
use BlockDist;
use Map;
use SysCTypes;

config const kProfilingThreads = true;


extern proc chpl_task_getId(): chpl_taskID_t;

require "util.h";
extern proc matvec_get_thread_id() : c_long;

inline proc getThreadId() { return matvec_get_thread_id():int; }

class MeasurementTable {
  var totalTime : [LocaleSpace dmapped Block(LocaleSpace)] map(int, real, true);
  var functionName : string;

  proc init(functionName : string) {
    this.functionName = functionName;
  }

  proc print() {
    writeln("Timings for ", functionName, ":");
    for loc in Locales {
      if (totalTime[loc.id].isEmpty()) { continue; }
      if (kProfilingThreads) {
        writeln("  locale ", loc.id, ": ", totalTime[loc.id]);
      }
      else {
        const total = max reduce totalTime[loc.id].values();
        writeln("  locale ", loc.id, ": ", total);
      }
    }
  }
}

record ScopedTimer {
  var table : borrowed MeasurementTable;
  var timer : Timer;

  inline proc init(ref table : MeasurementTable) {
    this.table = table.borrow();
    this.timer = new Timer();
    complete();
    timer.start();
  }

  inline proc stop() { timer.stop(); }

  inline proc deinit() {
    if (timer.running) { timer.stop(); }
    const myId = getThreadId(); // chpl_task_getId();
    const elapsed = timer.elapsed();
    ref myTable = table.totalTime[here.id];

    if (!myTable.add(myId, elapsed)) {
      record Updater {
        var elapsed : real;
        inline proc this(k, ref v) {
          v += elapsed;
          return none;
        }
      }
      myTable.update(myId, new Updater(elapsed));
    }
  }
}

inline proc getTimerFor(ref table : MeasurementTable) {
  return new ScopedTimer(table);
}
