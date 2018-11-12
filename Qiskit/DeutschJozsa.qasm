// Name of Experiment: DeutschJozsa v8

OPENQASM 2.0;
include "qelib1.inc";

qreg q[5];
creg c[5];

h q[0];
h q[1];
x q[2];
h q[2];
barrier q[0],q[1],q[2];
h q[0];
h q[1];
h q[2];
cx q[2],q[0];
cx q[2],q[1];
h q[0];
h q[1];
h q[2];
barrier q[0],q[1],q[2];
h q[2];
h q[0];
h q[1];
x q[2];
barrier q[0],q[1],q[2];
measure q[0] -> c[0];
measure q[1] -> c[1];
