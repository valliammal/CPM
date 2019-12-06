# CPM
This is for the critical Path. This work is based on the IEEE paper
We describe a new exact procedure for the discrete time=cost trade-off problem in deterministic 
activity-on-the-arc networks of the CPM type, where the duration of each activity is a discrete, 
nonincreasing function of the amount of a singleresource(money)committedtoit.
Theobjectiveistoconstructthecompleteandef®cienttime=costpro®leovertheset of feasible project durations. 
The procedure uses a horizon-varying approach based on the iterative optimal solution of the problem 
of minimising the sum of the resource use over all activities subject to the activity precedence 
constraints and a project deadline. This optimal solution is derived using a branch-and-bound procedure 
which computes lower bounds by making convex piecewise linear underestimations of the discrete time=cost 
trade-off curves of the activities to be used as input for an adapted version of the Fulkerson 
labelling algorithm for the linear time=cost trade-off problem. Branching involves the selection 
of an activity in order to partition its set of execution modes into two subsets which are used to 
derive improved convex piecewise linear underestimations. The procedure has been programmed in C++ under 
Windows NT and has been validated using a factorial experiment on a large set of randomly generated problem instances.

Please add the Boost libraries in the developed code
