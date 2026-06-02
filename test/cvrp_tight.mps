NAME          cvrp_tight
ROWS
 N  OBJ
 E  OUT1
 E  OUT2
 E  OUT3
 E  OUT4
 E  IN1
 E  IN2
 E  IN3
 E  IN4
 E  FLOW1
 E  FLOW2
 E  FLOW3
 E  FLOW4
 L  FLEET
 L  CAP0
 L  CAP1
 L  CAP2
 L  CAP3
 L  CAP4
 L  CAP5
 L  CAP6
 L  CAP7
 L  CAP8
 L  CAP9
 L  CAP10
 L  CAP11
 L  CAP12
 L  CAP13
 L  CAP14
 L  CAP15
 G  DLBO1
 L  DUBO1
 G  DLBO2
 L  DUBO2
 G  DLBO3
 L  DUBO3
 G  DLBO4
 L  DUBO4
COLUMNS
    MARK0000  'MARKER'                 'INTORG'
    x0_1          OBJ       15
    x0_1          IN1      1
    x0_1          FLOW1   -1
    x0_1          FLEET      1
    x0_1          CAP0    10
    x0_2          OBJ       20
    x0_2          IN2      1
    x0_2          FLOW2   -1
    x0_2          FLEET      1
    x0_2          CAP1    10
    x0_3          OBJ       25
    x0_3          IN3      1
    x0_3          FLOW3   -1
    x0_3          FLEET      1
    x0_3          CAP2    10
    x0_4          OBJ       30
    x0_4          IN4      1
    x0_4          FLOW4   -1
    x0_4          FLEET      1
    x0_4          CAP3    10
    x1_0          OBJ       15
    x1_0          OUT1     1
    x1_0          FLOW1    1
    x1_2          OBJ       15
    x1_2          OUT1     1
    x1_2          IN2      1
    x1_2          FLOW1    1
    x1_2          FLOW2   -1
    x1_2          CAP4    10
    x1_3          OBJ       20
    x1_3          OUT1     1
    x1_3          IN3      1
    x1_3          FLOW1    1
    x1_3          FLOW3   -1
    x1_3          CAP5    10
    x1_4          OBJ       25
    x1_4          OUT1     1
    x1_4          IN4      1
    x1_4          FLOW1    1
    x1_4          FLOW4   -1
    x1_4          CAP6    10
    x2_0          OBJ       20
    x2_0          OUT2     1
    x2_0          FLOW2    1
    x2_1          OBJ       15
    x2_1          OUT2     1
    x2_1          IN1      1
    x2_1          FLOW2    1
    x2_1          FLOW1   -1
    x2_1          CAP7    10
    x2_3          OBJ       15
    x2_3          OUT2     1
    x2_3          IN3      1
    x2_3          FLOW2    1
    x2_3          FLOW3   -1
    x2_3          CAP8    10
    x2_4          OBJ       20
    x2_4          OUT2     1
    x2_4          IN4      1
    x2_4          FLOW2    1
    x2_4          FLOW4   -1
    x2_4          CAP9    10
    x3_0          OBJ       25
    x3_0          OUT3     1
    x3_0          FLOW3    1
    x3_1          OBJ       20
    x3_1          OUT3     1
    x3_1          IN1      1
    x3_1          FLOW3    1
    x3_1          FLOW1   -1
    x3_1          CAP10    10
    x3_2          OBJ       15
    x3_2          OUT3     1
    x3_2          IN2      1
    x3_2          FLOW3    1
    x3_2          FLOW2   -1
    x3_2          CAP11    10
    x3_4          OBJ       15
    x3_4          OUT3     1
    x3_4          IN4      1
    x3_4          FLOW3    1
    x3_4          FLOW4   -1
    x3_4          CAP12    10
    x4_0          OBJ       30
    x4_0          OUT4     1
    x4_0          FLOW4    1
    x4_1          OBJ       25
    x4_1          OUT4     1
    x4_1          IN1      1
    x4_1          FLOW4    1
    x4_1          FLOW1   -1
    x4_1          CAP13    10
    x4_2          OBJ       20
    x4_2          OUT4     1
    x4_2          IN2      1
    x4_2          FLOW4    1
    x4_2          FLOW2   -1
    x4_2          CAP14    10
    x4_3          OBJ       15
    x4_3          OUT4     1
    x4_3          IN3      1
    x4_3          FLOW4    1
    x4_3          FLOW3   -1
    x4_3          CAP15    10
    MARK0000  'MARKER'                 'INTEND'
    u1            CAP0    -1
    u1            CAP4     1
    u1            CAP5     1
    u1            CAP6     1
    u1            CAP7    -1
    u1            CAP10    -1
    u1            CAP13    -1
    u1            DLBO1    1
    u1            DUBO1    1
    u2            CAP1    -1
    u2            CAP4    -1
    u2            CAP7     1
    u2            CAP8     1
    u2            CAP9     1
    u2            CAP11    -1
    u2            CAP14    -1
    u2            DLBO2    1
    u2            DUBO2    1
    u3            CAP2    -1
    u3            CAP5    -1
    u3            CAP8    -1
    u3            CAP10     1
    u3            CAP11     1
    u3            CAP12     1
    u3            CAP15    -1
    u3            DLBO3    1
    u3            DUBO3    1
    u4            CAP3    -1
    u4            CAP6    -1
    u4            CAP9    -1
    u4            CAP12    -1
    u4            CAP13     1
    u4            CAP14     1
    u4            CAP15     1
    u4            DLBO4    1
    u4            DUBO4    1
RHS
    RHS1      OUT1      1
    RHS1      IN1       1
    RHS1      OUT2      1
    RHS1      IN2       1
    RHS1      OUT3      1
    RHS1      IN3       1
    RHS1      OUT4      1
    RHS1      IN4       1
    RHS1      FLEET      2
    RHS1      CAP0     7
    RHS1      CAP1     6
    RHS1      CAP2     5
    RHS1      CAP3     4
    RHS1      CAP4     6
    RHS1      CAP5     5
    RHS1      CAP6     4
    RHS1      CAP7     7
    RHS1      CAP8     5
    RHS1      CAP9     4
    RHS1      CAP10     7
    RHS1      CAP11     6
    RHS1      CAP12     4
    RHS1      CAP13     7
    RHS1      CAP14     6
    RHS1      CAP15     5
    RHS1      DLBO1    3
    RHS1      DUBO1    10
    RHS1      DLBO2    4
    RHS1      DUBO2    10
    RHS1      DLBO3    5
    RHS1      DUBO3    10
    RHS1      DLBO4    6
    RHS1      DUBO4    10
BOUNDS
ENDATA
