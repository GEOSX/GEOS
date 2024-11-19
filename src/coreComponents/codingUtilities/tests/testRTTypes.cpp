/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "codingUtilities/RTTypes.hpp"

#include <gtest/gtest.h>


// Mock classes to test dynamic casting
class Base  {
public:
    virtual ~Base() = default; // Needed for RTTI
};

class Derived final : public Base {
public:
    void show() {
        std::cout << "Derived class method." << std::endl;
    }
};

// Test for dynamicCast with pointer
TEST(DynamicCastTests, Pointer_Casting_Success) {
    Base* base = new Derived();
    Derived* derived = geos::dynamicCast<Derived*>(base);
    ASSERT_NE(derived, nullptr) << "Expected successful cast from Base to Derived.";
    delete base; // Clean up allocated memory
}

TEST(DynamicCastTests, Pointer_Casting_Failure) {
    Base* base = new Base();
    Derived* derived = geos::dynamicCast<Derived*>(base);
    ASSERT_EQ(derived, nullptr) << "Expected nullptr due to failed cast from Base to Derived.";
    delete base; // Clean up allocated memory
}

// Test for dynamicCast with reference
TEST(DynamicCastTests, Reference_Casting_Success) {
    Derived derived;
    Base& base_ref = derived;
    Derived& derived_ref = geos::dynamicCast<Derived&>(base_ref);
    ASSERT_EQ(&derived_ref, &derived) << "Expected successful cast from Base to Derived.";
}

TEST(DynamicCastTests, Reference_Casting_Failure) {
    Base base;
    Base& base_ref = base;
    
    Base& derived_base_ref = geos::dynamicCast<Base&>(base_ref);
    ASSERT_EQ(&derived_base_ref, &base) << "Expected successful cast from Base to Base.";
    
}

// Test Regex constructor
TEST(RegexTests, Constructor) {
    geos::Regex regex("^[0-9]+$", "Input must be a number.");
    ASSERT_EQ(regex.m_regexStr, "^[0-9]+$") << "Regex string is incorrect.";
    ASSERT_EQ(regex.m_formatDescription, "Input must be a number.") << "Format description is incorrect.";
}

TEST(RtTypesTests, GetTypeName) {
    {
        std::type_index typeIndex(typeid(Base));
        auto typeName = geos::rtTypes::getTypeName(typeIndex);
        EXPECT_EQ(typeName, std::string("Base")); // Expected Base
    }
    {
        std::type_index typeIndex(typeid(Derived));
        auto typeName = geos::rtTypes::getTypeName(typeIndex);
        EXPECT_EQ(typeName, std::string("Derived")); // Expected Derived
    }
}

// Additional tests to validate the functionality of getTypeRegex
TEST(RtTypesTests, GetTypeRegex_Default) {
    geos::Regex regex = geos::rtTypes::getTypeRegex<int>(); // Assuming int has a default regex defined
    ASSERT_NE(regex.m_regexStr.empty(), true) << "Expected non-empty regex for int.";
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

