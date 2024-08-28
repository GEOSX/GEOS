/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file testReverseCutHillMcKeeOrdering.cpp
 */

#include "codingUtilities/UnitTestUtilities.hpp"
#include "linearAlgebra/utilities/ReverseCutHillMcKeeOrdering.hpp"

#include <gtest/gtest.h>

using namespace geos;

TEST( ReverseCutHillMcKeeOrderingTest, reorder )
{
  // the data in this test comes from PoroElastic_staircase_co2_3d.xml
  // it is extracted of the flow dof numbers (MPI rank 1 in a 2-rank simulation)

  localIndex constexpr numRows = 144;
  localIndex constexpr rankOffset = 144;
  localIndex constexpr numNonZeros = 876;

  array1d< localIndex > permutation( numRows );
  localIndex const expectedPermutation[numRows]
  { 16, 10, 17, 14, 4, 8, 11, 34, 15, 12, 2, 5, 6, 9, 28, 35, 32, 13,
    88, 0, 3, 22, 82, 7, 26, 29, 106, 33, 30, 89, 86, 76, 1, 20, 23, 80, 83,
    24, 27, 100, 107, 104, 31, 52, 87, 84, 74, 77, 18, 21, 94, 78, 81, 46,
    25, 98, 101, 124, 105, 102, 53, 50, 85, 72, 75, 40, 19, 92, 95, 79, 44,
    47, 96, 99, 118, 125, 122, 103, 70, 51, 48, 73, 38, 41, 90, 93, 112, 42,
    45, 64, 97, 116, 119, 123, 120, 71, 68, 49, 36, 39, 58, 91, 110, 113, 43,
    62, 65, 114, 117, 121, 142, 69, 66, 37, 56, 59, 108, 111, 60, 63, 136,
    115, 143, 140, 67, 54, 57, 130, 109, 61, 134, 137, 141, 138, 55, 128, 131,
    132, 135, 139, 126, 129, 133, 127 };
  localIndex const offsets[numRows+1]
  { 0, 6, 13, 19, 26, 31, 37, 43, 50, 56, 63, 68, 74, 79, 85, 90, 96,
    100, 105, 112, 119, 126, 133, 139, 145, 152, 159, 166, 173, 179, 185, 191,
    197, 203, 209, 214, 219, 225, 231, 238, 245, 252, 259, 265, 271, 278, 285,
    292, 299, 304, 309, 315, 321, 327, 333, 339, 345, 352, 359, 366, 373, 379,
    385, 392, 399, 406, 413, 418, 423, 429, 435, 441, 447, 452, 458, 464, 471,
    477, 484, 489, 495, 501, 508, 514, 521, 525, 530, 535, 541, 546, 552, 559,
    566, 573, 580, 586, 592, 599, 606, 613, 620, 626, 632, 638, 644, 650, 656,
    661, 666, 673, 679, 686, 692, 698, 703, 710, 716, 723, 729, 735, 740, 746,
    751, 757, 762, 767, 771, 777, 782, 789, 795, 802, 808, 814, 819, 826, 832,
    839, 845, 850, 854, 860, 865, 871, 876  };
  globalIndex const columns[numNonZeros+1]
  { 12, 144, 145, 146, 150, 220, 13, 144, 145, 147, 151, 162, 221, 14, 144, 146,
    147, 148, 152, 15, 145, 146, 147, 149, 153, 164, 16, 146, 148, 149, 154, 17,
    147, 148, 149, 155, 166, 144, 150, 151, 152, 156, 226, 145, 150, 151, 153, 157,
    168, 227, 146, 150, 152, 153, 154, 158, 147, 151, 152, 153, 155, 159, 170, 148,
    152, 154, 155, 160, 149, 153, 154, 155, 161, 172, 150, 156, 157, 158, 232, 151,
    156, 157, 159, 174, 233, 152, 156, 158, 159, 160, 153, 157, 158, 159, 161, 176,
    154, 158, 160, 161, 155, 159, 160, 161, 178, 120, 145, 162, 163, 164, 168, 184,
    121, 162, 163, 165, 169, 185, 234, 122, 147, 162, 164, 165, 166, 170, 123, 163,
    164, 165, 167, 171, 236, 124, 149, 164, 166, 167, 172, 125, 165, 166, 167, 173,
    238, 151, 162, 168, 169, 170, 174, 190, 163, 168, 169, 171, 175, 191, 240, 153,
    164, 168, 170, 171, 172, 176, 165, 169, 170, 171, 173, 177, 242, 155, 166, 170,
    172, 173, 178, 167, 171, 172, 173, 179, 244, 157, 168, 174, 175, 176, 196, 169,
    174, 175, 177, 197, 246, 159, 170, 174, 176, 177, 178, 171, 175, 176, 177, 179,
    248, 161, 172, 176, 178, 179, 173, 177, 178, 179, 250, 102, 180, 181, 182, 186,
    217, 103, 180, 181, 183, 187, 198, 104, 180, 182, 183, 184, 188, 219, 105, 181,
    182, 183, 185, 189, 200, 106, 162, 182, 184, 185, 190, 221, 107, 163, 183, 184,
    185, 191, 202, 180, 186, 187, 188, 192, 223, 181, 186, 187, 189, 193, 204, 182,
    186, 188, 189, 190, 194, 225, 183, 187, 188, 189, 191, 195, 206, 168, 184, 188,
    190, 191, 196, 227, 169, 185, 189, 190, 191, 197, 208, 186, 192, 193, 194, 229,
    187, 192, 193, 195, 210, 188, 192, 194, 195, 196, 231, 189, 193, 194, 195, 197,
    212, 174, 190, 194, 196, 197, 233, 175, 191, 195, 196, 197, 214, 30, 181, 198,
    199, 200, 204, 31, 198, 199, 201, 205, 270, 32, 183, 198, 200, 201, 202, 206,
    33, 199, 200, 201, 203, 207, 272, 34, 185, 200, 202, 203, 208, 234, 35, 201,
    202, 203, 209, 235, 274, 187, 198, 204, 205, 206, 210, 199, 204, 205, 207, 211,
    276, 189, 200, 204, 206, 207, 208, 212, 201, 205, 206, 207, 209, 213, 278, 191,
    202, 206, 208, 209, 214, 240, 203, 207, 208, 209, 215, 241, 280, 193, 204, 210,
    211, 212, 205, 210, 211, 213, 282, 195, 206, 210, 212, 213, 214, 207, 211, 212,
    213, 215, 284, 197, 208, 212, 214, 215, 246, 209, 213, 214, 215, 247, 286, 84,
    216, 217, 218, 222, 85, 180, 216, 217, 219, 223, 86, 216, 218, 219, 220, 224,
    87, 182, 217, 218, 219, 221, 225, 88, 144, 218, 220, 221, 226, 89, 145, 184,
    219, 220, 221, 227, 216, 222, 223, 224, 228, 186, 217, 222, 223, 225, 229, 218,
    222, 224, 225, 226, 230, 188, 219, 223, 224, 225, 227, 231, 150, 220, 224, 226,
    227, 232, 151, 190, 221, 225, 226, 227, 233, 222, 228, 229, 230, 192, 223, 228,
    229, 231, 224, 228, 230, 231, 232, 194, 225, 229, 230, 231, 233, 156, 226, 230,
    232, 233, 157, 196, 227, 231, 232, 233, 138, 163, 202, 234, 235, 236, 240, 139,
    203, 234, 235, 237, 241, 252, 140, 165, 234, 236, 237, 238, 242, 141, 235, 236,
    237, 239, 243, 254, 142, 167, 236, 238, 239, 244, 143, 237, 238, 239, 245, 256,
    169, 208, 234, 240, 241, 242, 246, 209, 235, 240, 241, 243, 247, 258, 171, 236,
    240, 242, 243, 244, 248, 237, 241, 242, 243, 245, 249, 260, 173, 238, 242, 244,
    245, 250, 239, 243, 244, 245, 251, 262, 175, 214, 240, 246, 247, 248, 215, 241,
    246, 247, 249, 264, 177, 242, 246, 248, 249, 250, 243, 247, 248, 249, 251, 266,
    179, 244, 248, 250, 251, 245, 249, 250, 251, 268, 66, 235, 252, 253, 254, 258,
    274, 67, 252, 253, 255, 259, 275, 68, 237, 252, 254, 255, 256, 260, 69, 253, 254,
    255, 257, 261, 70, 239, 254, 256, 257, 262, 71, 255, 256, 257, 263, 241, 252, 258,
    259, 260, 264, 280, 253, 258, 259, 261, 265, 281, 243, 254, 258, 260, 261, 262,
    266, 255, 259, 260, 261, 263, 267, 245, 256, 260, 262, 263, 268, 257, 261, 262,
    263, 269, 247, 258, 264, 265, 266, 286, 259, 264, 265, 267, 287, 249, 260, 264,
    266, 267, 268, 261, 265, 266, 267, 269, 251, 262, 266, 268, 269, 263, 267, 268,
    269, 48, 199, 270, 271, 272, 276, 49, 270, 271, 273, 277, 50, 201, 270, 272, 273,
    274, 278, 51, 271, 272, 273, 275, 279, 52, 203, 252, 272, 274, 275, 280, 53, 253,
    273, 274, 275, 281, 205, 270, 276, 277, 278, 282, 271, 276, 277, 279, 283, 207,
    272, 276, 278, 279, 280, 284, 273, 277, 278, 279, 281, 285, 209, 258, 274, 278,
    280, 281, 286, 259, 275, 279, 280, 281, 287, 211, 276, 282, 283, 284, 277, 282,
    283, 285, 213, 278, 282, 284, 285, 286, 279, 283, 284, 285, 287, 215, 264, 280,
    284, 286, 287, 265, 281, 285, 286, 287 };

  reverseCutHillMcKeeOrdering::
    computePermutation( offsets, columns, rankOffset, permutation );

  for( localIndex i = 0; i < numRows; ++i )
  {
    EXPECT_EQ( permutation[i], expectedPermutation[i] );
  }
}
