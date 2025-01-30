//
//  DerivedClassFinal.hpp
//  testRTTypes
//
//  Created by Omar Duran on 12/16/24.
//

#ifndef DerivedClassFinal_hpp
#define DerivedClassFinal_hpp

#include <stdio.h>

class Derived final : public Base
{
public:
    explicit Derived();
    
    virtual ~Derived() noexcept override {
        
    };
};

#endif /* DerivedClassFinal_hpp */
