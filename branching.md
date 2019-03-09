# Branching

CBC's concept of branching is based on the idea of an "object".
An object has (i) a feasible region, (ii) can be evaluated for infeasibility, (iii) can be branched on, e.g., a method of generating a branching object, which defines an up branch and a down branch, and (iv) allows comparsion of the effect of branching.
Instances of objects include:
- `CbcSimpleInteger`,
- `CbcSimpleIntegerPseudoCosts`,
- `CbcClique`,
- `CbcSOS` (type 1 and 2),
- `CbcFollowOn`, and
- `CbcLotsize`.

## Pseudo Cost Branching

If the user declares variables as integer but does no more, then CBC
will treat them as simple integer variables. In many cases the user
would like to do some more fine tuning. This section shows how to create
integer variables with pseudo costs. When pseudo costs are given then it
is assumed that if a variable is at 1.3 then the cost of branching that
variable down will be 0.3 times the down pseudo cost and the cost of
branching up would be 0.7 times the up pseudo cost. Pseudo costs can be
used both for branching and for choosing a node. The full code is in
`longthin.cpp` located in the CBC Samples directory.

The idea is simple for set covering problems. Branching up gets us much
closer to an integer solution so we will encourage that direction by
branch up if variable value > 0.333333. The expected cost of going up
obviously depends on the cost of the variable. The pseudo costs are
choosen to reflect that fact.

```
int iColumn;
int numberColumns = solver3->getNumCols();
// do pseudo costs
CbcObject ** objects = new CbcObject * [numberColumns];
// Point to objective
const double * objective = model.getObjCoefficients();
int numberIntegers=0;
for (iColumn=0;iColumn<numberColumns;iColumn++) {
  if (solver3->isInteger(iColumn)) {
    double cost = objective[iColumn];
    CbcSimpleIntegerPseudoCost * newObject =
      new CbcSimpleIntegerPseudoCost(&model,numberIntegers,iColumn,
                                     2.0*cost,cost);
    newObject->setMethod(3);
    objects[numberIntegers++]= newObject;
  }
}
// Now add in objects (they will replace simple integers)
model.addObjects(numberIntegers,objects);
for (iColumn=0;iColumn<numberIntegers;iColumn++)
  delete objects[iColumn];
delete [] objects;
```

This code also tries to give more importance
to variables with more coefficients. Whether this sort of thing is
worthwhile should be the subject of experimentation.

## Follow-On Branching

In crew scheduling, the problems are long and thin. A problem may have a
few rows but many thousands of variables. Branching a variable to 1 is
very powerful as it fixes many other variables to zero, but branching to
zero is very weak as thousands of variables can increase from zero. In
crew scheduling problems, each constraint is a flight leg, e.g., JFK
airport to DFW airport. From DFW there may be several flights the crew
could take next - suppose one flight is the 9:30 flight from DFW to LAX
airport. A binary branch is that the crew arriving at DFW either take
the 9:30 flight to LAX or they don't. This "follow-on" branching does
not fix individual variables. Instead this branching divides all the
variables with entries in the JFK-DFW constraint into two groups - those
with entries in the DFW-LAX constraint and those without entries.

The full sample code for follow-on brancing is in `crew.cpp` located in
the CBC Samples directory. In this case, the
simple integer variables are left which may be necessary if other sorts
of constraints exist. Follow-on branching rules are to be considered
first, so the priorities are set to indicated the follow-on rules take
precedence. Priority 1 is the highest priority.

```
int iColumn;
int numberColumns = solver3->getNumCols();
/* We are going to add a single follow-on object but we
   want to give low priority to existing integers
   As the default priority is 1000 we don't actually need to give
   integer priorities but it is here to show how.
*/
// Normal integer priorities
int * priority = new int [numberColumns];
int numberIntegers=0;
for (iColumn=0;iColumn<numberColumns;iColumn++) {
  if (solver3->isInteger(iColumn)) {
    priority[numberIntegers++]= 100; // low priority
  }
}
/* Second parameter is set to true for objects,
   and false for integers. This indicates integers */
model.passInPriorities(priority,false);
delete [] priority;
/* Add in objects before we can give them a priority.
   In this case just one object
   - but it shows the general method
*/
CbcObject ** objects = new CbcObject * [1];
objects[0]=new CbcFollowOn(&model);
model.addObjects(1,objects);
delete objects[0];
delete [] objects;
// High priority
int followPriority=1;
model.passInPriorities(&followPriority,true);
```
