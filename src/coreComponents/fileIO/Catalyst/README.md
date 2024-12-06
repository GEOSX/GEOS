## Reason for this module

When trying to input the `conduit::Node` used by the simulation into the `conduit_cpp::Node` needed by `catalyst` we run into an incompatibility between the two conduit libraries which cannot coexist in the same translation unit. As of the time of writing this README, `catalyst` does not have any ready-made solution for this type of collision with simulation codes.

The solution here is to create a library in this module with the following specifications:
- Provides an interface class defining an API for interacting with generic `conduit` nodes
- Provides utilities to convert any implementation of the above interface into a `conduit_cpp::Node` for passing to `catalyst`
- Links to the `conduit` provided by `catalyst` in its own translation unit
- Provides its own `initialize`, `execute` and `finalize` methods taking care of the conversion from the simulation code's `conduit` and `catalyst`'s version

As such, for any simulation code already instrumented with `conduit` blueprints, this design pattern enables coexistance and conversion between the two `conduit` versions. The simulation code has but to provide an implementation of the interface in this module (in a separate library that links to this one) and pass it to the specific initialization, execution and finalization methods implemented in this library. The interface node acts as a buffer layer between the two `conduit` libraries and separates the translation units cleanly.
