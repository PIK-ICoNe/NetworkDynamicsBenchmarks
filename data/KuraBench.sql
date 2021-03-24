CREATE TABLE edges (system integer NOT NULL, id integer NOT NULL, source integer, dest integer, PRIMARY KEY (system, id));
CREATE TABLE IF NOT EXISTS "systems" (
	"id"	INTEGER,
	"name"	TEXT,
	"nodes"	INTEGER,
	"edges"	INTEGER,
	"q"	NUMERIC,
	"c"	NUMERIC,
	"coupling_constant"	NUMERIC,
	PRIMARY KEY("id")
);
CREATE TABLE IF NOT EXISTS "nodes" (
	"system"	integer NOT NULL,
	"id"	integer NOT NULL,
	"omega"	DOUBLE,
	"initialCondition"	DOUBLE,
	FOREIGN KEY("system") REFERENCES "systems"("id"),
	PRIMARY KEY("system","id")
);
CREATE TABLE IF NOT EXISTS "states" (
	"experimentid"	INTEGER NOT NULL,
	"time"	NUMERIC,
	"idx"	INTEGER,
	"node"	INTEGER,
	"state"	DOUBLE,
	PRIMARY KEY("experimentid","idx","time"),
	FOREIGN KEY("experimentid") REFERENCES "experiments"("experimentid")
);
CREATE TABLE IF NOT EXISTS "experiments" (
	"experimentid"	INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
	"runtime"	TEXT,
	"datetime"	TEXT,
	"githash"	TEXT,
	"gitchanges"	TEXT,
	"commandline"	TEXT,
	"host"	TEXT,
	"system"	INTEGER,
	"solver"	TEXT,
	"atol"	NUMERIC,
	"rtol"	NUMERIC,
	"tstart"	INTEGER,
	"tend"	INTEGER,
	"nfcn"	INTEGER,
	"njac"	INTEGER,
	"nstep"	INTEGER,
	"naccept"	INTEGER,
	"nreject"	INTEGER,
	"ndec"	INTEGER,
	"nsol"	INTEGER,
	"tstartup"	NUMERIC,
	"tload"	NUMERIC,
	"tjit"	NUMERIC,
	"tintegration"	NUMERIC,
	"ttotal"	NUMERIC,
	"peak_mem"	NUMERIC,
	"err_v_ref"	NUMERIC,
	"refeid"	INTEGER,
	"notes"	TEXT,
	FOREIGN KEY("system") REFERENCES "systems"("id")
);
