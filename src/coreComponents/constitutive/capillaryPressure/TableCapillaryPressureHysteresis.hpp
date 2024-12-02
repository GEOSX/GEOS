/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


#ifndef GEOS_CONSTITUTIVE_TABLECAPILLARYPRESSUREHYSTERESIS_HPP
#define GEOS_CONSTITUTIVE_TABLECAPILLARYPRESSUREHYSTERESIS_HPP

#include "constitutive/capillaryPressure/CapillaryPressureBase.hpp"
#include "functions/TableFunction.hpp"

#include "constitutive/KilloughHysteresis.hpp"
#include "CapillaryPressureFields.hpp"

namespace geos {


    namespace constitutive {

        class TableCapillaryPressureHysteresis : public CapillaryPressureBase {
//  /// useful constant
//  static constexpr real64 CAP_INF = 1e9;
////          std::numeric_limits< real64 >::max();
//  static constexpr real64 CAP_INF_DERIV = 1e9;
////          std::numeric_limits< real64 >::max();

            typedef fields::cappres::ModeIndexType ModeIndexType;

        public:

            /// order of the phase properties for three-phase flow
            struct ThreePhasePairPhaseType {
                enum : integer {
                    INTERMEDIATE_WETTING = 0,   ///< index for intermediate-wetting
                    INTERMEDIATE_NONWETTING = 1 ///< index for intermediate-non-wetting
                };
            };


            TableCapillaryPressureHysteresis(std::string const &name,
                                             dataRepository::Group *const parent);

            static std::string catalogName() { return "TableCapillaryPressureHysteresis"; }

            virtual string getCatalogName() const override { return catalogName(); }

            ///Kernel
            class KernelWrapper final : public CapillaryPressureBaseUpdate {
            public:

                KernelWrapper(
                        arrayView1d<TableFunction::KernelWrapper const> const &wettingNonWettingCapillaryPressureKernelWrappers,
                        arrayView1d<TableFunction::KernelWrapper const> const &wettingIntermediateCapillaryPressureKernelWrappers,
                        arrayView1d<TableFunction::KernelWrapper const> const &nonWettingIntermediateCapillaryPressureKernelWrappers,
                        arrayView1d<integer const> const &phaseHasHysteresis,
                        arrayView1d<real64 const> const &landParam,
                        real64 const &jerauldParam_a,
                        real64 const &jerauldParam_b,
                        real64 const &killoughCurvaturePcParameter,
                        real64 const &phaseIntermediateMinVolFraction,
                        KilloughHysteresis::HysteresisCurve const &wettingCurve,
                        KilloughHysteresis::HysteresisCurve const &nonWettingCurve,
                        arrayView2d<real64 const, compflow::USD_PHASE> const &phaseMinHistoricalVolFraction,
                        arrayView2d<real64 const, compflow::USD_PHASE> const &phaseMaxHistoricalVolFraction,
                        arrayView1d<integer const> const &phaseTypes,
                        arrayView1d<integer const> const &phaseOrder,
                        arrayView1d<integer> const &mode,
                        arrayView3d<real64, relperm::USD_RELPERM> const &phaseTrappedVolFrac,
                        arrayView3d<real64, relperm::USD_RELPERM> const &phaseCapPressure,
                        arrayView4d<real64, relperm::USD_RELPERM_DS> const &dPhaseCapPressure_dPhaseVolFrac);

                //actual workers
                GEOS_HOST_DEVICE
                void computeBoundCapillaryPressure(TableFunction::KernelWrapper const &drainageCapPressureWrapper,
                                                   real64 const &phaseVolFraction,
                                                   real64 &phaseCapPressure,
                                                   real64 &dPhaseCapPressure_dPhaseVolFrac) const;

                GEOS_HOST_DEVICE
                void
                computeImbibitionWettingCapillaryPressure(
                        const arrayView1d<const TableFunction::KernelWrapper> &wettingKernelWapper,
                        const KilloughHysteresis::HysteresisCurve &wettingCurve,
                        const KilloughHysteresis::HysteresisCurve &nonWettingCurve,
                        const geos::real64 &landParam,
                        const geos::real64 &phaseVolFraction,
                        const geos::real64 &phaseMinHistoricalVolFraction,
                        geos::real64 &phaseTrappedVolFrac,
                        geos::real64 &phaseCapPressure,
                        geos::real64 &dPhaseCapPressure_dPhaseVolFrac,
                        const ModeIndexType &mode) const;

                //two phase flow overload
                GEOS_HOST_DEVICE
                void
                computeImbibitionWettingCapillaryPressure(
                        const arrayView1d<const TableFunction::KernelWrapper> &wettingKernelWapper,
                        const KilloughHysteresis::HysteresisCurve &wettingCurve,
                        const geos::real64 &landParam,
                        const geos::real64 &phaseVolFraction,
                        const geos::real64 &phaseMinHistoricalVolFraction,
                        geos::real64 &phaseTrappedVolFrac,
                        geos::real64 &phaseCapPressure,
                        geos::real64 &dPhaseCapPressure_dPhaseVolFrac,
                        const ModeIndexType &mode) const;

                GEOS_HOST_DEVICE
                void
                computeImbibitionNonWettingCapillaryPressure(
                        const arrayView1d<const TableFunction::KernelWrapper> &nonWettingKernelWrapper,
                        const KilloughHysteresis::HysteresisCurve &nonWettingCurve,
                        const KilloughHysteresis::HysteresisCurve &wettingCurve,
                        const geos::real64 &landParam,
                        const geos::real64 &phaseVolFraction,
                        const geos::real64 &phaseMaxHistoricalVolFraction,
                        geos::real64 &phaseTrappedVolFrac,
                        geos::real64 &phaseCapPressure,
                        geos::real64 &dPhaseCapPressure_dPhaseVolFrac,
                        const ModeIndexType &mode) const;

                //2phase flow overload
                GEOS_HOST_DEVICE
                void
                computeImbibitionNonWettingCapillaryPressure(
                        const arrayView1d<const TableFunction::KernelWrapper> &nonWettingKernelWrapper,
                        const KilloughHysteresis::HysteresisCurve &nonWettingCurve,
                        const geos::real64 &landParam,
                        const geos::real64 &phaseVolFraction,
                        const geos::real64 &phaseMaxHistoricalVolFraction,
                        geos::real64 &phaseTrappedVolFrac,
                        geos::real64 &phaseCapPressure,
                        geos::real64 &dPhaseCapPressure_dPhaseVolFrac,
                        const ModeIndexType &mode) const;



                //wrapper call wrt number of phase
                GEOS_HOST_DEVICE
                void computeTwoPhaseWetting(integer const ipWetting,
                                            integer const ipNonWetting,
                                            arraySlice1d<real64 const, compflow::USD_PHASE - 1> const &phaseVolFraction,
                                            arraySlice1d<real64 const,
                                                    compflow::USD_PHASE - 1> const &phaseMaxHistoricalVolFraction,
                                            arraySlice1d<real64 const,
                                                    compflow::USD_PHASE - 1> const &phaseMinHistoricalVolFraction,
                                            arraySlice1d<real64, relperm::USD_RELPERM - 2> const &phaseTrappedVolFrac,
                                            arraySlice1d<real64, relperm::USD_RELPERM - 2> const &phaseCapPressure,
                                            arraySlice2d<real64,
                                                    relperm::USD_RELPERM_DS - 2> const &dPhaseCapPressure_dPhaseVolFrac,
                                            ModeIndexType &mode) const;


                GEOS_HOST_DEVICE
                void computeTwoPhaseNonWetting(integer const ipWetting,
                                               integer const ipNonWetting,
                                               arraySlice1d<real64 const,
                                                       compflow::USD_PHASE - 1> const &phaseVolFraction,
                                               arraySlice1d<real64 const,
                                                       compflow::USD_PHASE - 1> const &phaseMaxHistoricalVolFraction,
                                               arraySlice1d<real64 const,
                                                       compflow::USD_PHASE - 1> const &phaseMinHistoricalVolFraction,
                                               arraySlice1d<real64,
                                                       relperm::USD_RELPERM - 2> const &phaseTrappedVolFrac,
                                               arraySlice1d<real64, relperm::USD_RELPERM - 2> const &phaseCapPressure,
                                               arraySlice2d<real64, relperm::USD_RELPERM_DS -
                                                                    2> const &dPhaseCapPressure_dPhaseVolFrac,
                                               ModeIndexType &mode) const;

                GEOS_HOST_DEVICE
                void computeThreePhase(integer const ipWetting,
                                       integer const ipInter,
                                       integer const ipNonWetting,
                                       arraySlice1d<real64 const, compflow::USD_PHASE - 1> const &phaseVolFraction,
                                       arraySlice1d<real64 const,
                                               compflow::USD_PHASE - 1> const &phaseMaxHistoricalVolFraction,
                                       arraySlice1d<real64 const,
                                               compflow::USD_PHASE - 1> const &phaseMinHistoricalVolFraction,
                                       arraySlice1d<real64, relperm::USD_RELPERM - 2> const &phaseTrappedVolFrac,
                                       arraySlice1d<real64, relperm::USD_RELPERM - 2> const &phaseCapPressure,
                                       arraySlice2d<real64,
                                               relperm::USD_RELPERM_DS - 2> const &dPhaseCapPressure_dPhaseVolFrac,
                                       ModeIndexType &mode) const;

                //uppermost call-wrappers
                GEOS_HOST_DEVICE
                virtual void compute(arraySlice1d<real64 const, compflow::USD_PHASE - 1> const &phaseVolFraction,
                                     arraySlice1d<real64 const,
                                             compflow::USD_PHASE - 1> const &phaseMaxHistoricalVolFraction,
                                     arraySlice1d<real64 const,
                                             compflow::USD_PHASE - 1> const &phaseMinHistoricalVolFraction,
                                     arraySlice1d<real64, cappres::USD_CAPPRES - 2> const &phaseTrappedVolFrac,
                                     arraySlice1d<real64, cappres::USD_CAPPRES - 2> const &phaseCapPressure,
                                     arraySlice2d<real64,
                                             cappres::USD_CAPPRES_DS - 2> const &dPhaseCapPressure_dPhaseVolFrac,
                                     ModeIndexType &mode) const;

                GEOS_HOST_DEVICE
                virtual void update(localIndex const k,
                                    localIndex const q,
                                    arraySlice1d<real64 const,
                                            compflow::USD_PHASE - 1> const &phaseVolFraction) const override;

            private:


                static constexpr real64 flowReversalBuffer = KilloughHysteresis::flowReversalBuffer;
//    ModeIndexType& m_mode;

                //2p
                arrayView1d<TableFunction::KernelWrapper const> const m_wettingNonWettingCapillaryPressureKernelWrappers;
                //3p
                arrayView1d<TableFunction::KernelWrapper const> const m_wettingIntermediateCapillaryPressureKernelWrappers;
                arrayView1d<TableFunction::KernelWrapper const> const m_nonWettingIntermediateCapillaryPressureKernelWrappers;

                ///Land Coeff
                arrayView1d<integer const> m_phaseHasHysteresis;
                arrayView1d<real64 const> m_landParam;

                /// Parameter a introduced by Jerauld in the Land model
                const real64 m_jerauldParam_a;

                /// Parameter b introduced by Jerauld in the Land model
                const real64 m_jerauldParam_b;

                /// Curvature parameter in Killough wetting phase hysteresis (enpoints curvatures)
                const real64 m_killoughCurvatureParamCapPres;

                /// needed in 3p-wetting hysteresis as we need to get the max accessible pore space
                real64 const m_phaseIntermediateMinVolFraction;

                KilloughHysteresis::HysteresisCurve const m_wettingCurve;
                KilloughHysteresis::HysteresisCurve const m_nonWettingCurve;

                /// Minimum historical phase volume fraction for each phase
                arrayView2d<real64 const, compflow::USD_PHASE> m_phaseMinHistoricalVolFraction;

                /// Maximum historical phase volume fraction for each phase
                arrayView2d<real64 const, compflow::USD_PHASE> m_phaseMaxHistoricalVolFraction;

                // Drainage / Imbibition flags cellwise
                arrayView1d<ModeIndexType> m_mode;

            };

            /**
             * @brief Create an update kernel wrapper.
             * @return the wrapper
             */
            KernelWrapper createKernelWrapper();

            //might need it to be virtual one level higher --> from Killough/Hysteresis common class
            virtual void saveConvergedPhaseVolFractionState(
                    arrayView2d<real64 const, compflow::USD_PHASE> const &phaseVolFraction) const override;


            struct viewKeyStruct : CapillaryPressureBase::viewKeyStruct {


                ///Land Coeff
                static constexpr char const *landParameterString() { return "landParameter"; }

                ///flag
                static constexpr char const *phaseHasHysteresisString() { return "phaseHasHysteresis"; }

                ///and packed curves data struct
                static constexpr char const *wettingCurveString() { return "wettingCurve"; };

                static constexpr char const *nonWettingCurveString() { return "nonWettingCurve"; };


                ///tables and assoc. wrappers
                //2phase
                static constexpr char const *
                drainageWettingNonWettingCapPresTableNameString() { return "drainageWettingNonWettingCapPressureTableName"; }

                static constexpr char const *
                imbibitionWettingNonWettingCapPresTableNameString() { return "imbibitionWettingNonWettingCapPressureTableName"; }

                //3phase
                static constexpr char const *
                drainageWettingIntermediateCapPresTableNameString() { return "drainageWettingIntermediateCapPressureTableName"; }

                static constexpr char const *
                drainageNonWettingIntermediateCapPresTableNameString() { return "drainageNonWettingIntermediateCapPressureTableName"; }

                static constexpr char const *
                imbibitionWettingIntermediateCapPresTableNameString() { return "imbibitionWettingIntermediateCapPressureTableName"; }

                static constexpr char const *
                imbibitionNonWettingIntermediateCapPresTableNameString() { return "imbibitionNonWettingIntermediateCapPressureTableName"; }

                static constexpr char const *
                wettingNonWettingCapillaryPressureKernelWrappersString() { return "wettingNonWettingCapillaryPressureKernelWrappers"; }

                static constexpr char const *
                wettingIntermediateCapillaryPressureKernelWrappersString() { return "wettingIntermediateCapillaryPressureKernelWrappers"; }

                static constexpr char const *
                nonWettingIntermediateCapillaryPressureKernelWrappersString() { return "nonWettingIntermediateCapillaryPressureKernelWrappers"; }

                //misc
                static constexpr char const *
                phaseIntermediateMinVolFractionString() { return "phaseIntermediateMinVolFraction"; }
                //to decide wheter drainage/drainage to imbibition or imbibition/imbibition to drainage
            };


        private:
            virtual void postProcessInput() override;

            virtual void initializePreSubGroups() override;

            void resizeFields(localIndex const size,
                              localIndex const numPts) override;


            /**
             * @brief Create all the table kernel wrappers needed for the simulation (for all the phases present)
             */
            void createAllTableKernelWrappers();

            /**
             * @brief Compute the Land coefficient for the wetting and non-wetting phases
             */
            void computeLandCoefficient();

            ///data members


            //TODO impl
//  array1d< integer >  m_tCurveOption;

            KilloughHysteresis::HysteresisCurve m_wettingCurve;
            KilloughHysteresis::HysteresisCurve m_nonWettingCurve;

            ///tables
            //2p
            string m_drainageWettingNonWettingCapPresTableName;
            string m_imbibitionWettingNonWettingCapPresTableName;
            //3p
            string m_drainageWettingIntermediateCapPresTableName;
            string m_drainageNonWettingIntermediateCapPresTableName;
            string m_imbibitionWettingIntermediateCapPresTableName;
            string m_imbibitionNonWettingIntermediateCapPresTableName;
            // kernel wrappers
            /// Imbibition kernel wrappers for relative permeabilities in the following order:
            /// 0- drainage
            /// 1- imbibition (cf. struct ModeIndexType)
            //2p
            array1d<TableFunction::KernelWrapper> m_wettingNonWettingCapillaryPressureKernelWrappers;
            //3p
            array1d<TableFunction::KernelWrapper> m_wettingIntermediateCapillaryPressureKernelWrappers;
            array1d<TableFunction::KernelWrapper> m_nonWettingIntermediateCapillaryPressureKernelWrappers;


            /// Flag to specify whether the phase has hysteresis or not (deduced from table input)
            array1d<integer> m_phaseHasHysteresis;

            /// Trapping parameter from the Land model (typically called C)
            array1d<real64> m_landParam;

            /// Parameter a introduced by Jerauld in the Land model
            real64 m_jerauldParam_a;

            /// Parameter b introduced by Jerauld in the Land model
            real64 m_jerauldParam_b;

            /// Curvature parameter in Killough wetting phase hysteresis (Scanning curves curvatures)
            real64 m_killoughCurvatureParamCapPres;

            /// Cell-wise status imbibition, imbibitioon_to_drainage, ... etc
            array1d<integer> m_mode;

            // Max historical saturations
            /// Minimum historical phase volume fraction for each phase
            array2d<real64, compflow::LAYOUT_PHASE> m_phaseMinHistoricalVolFraction;

            /// Maximum historical phase volume fraction for each phase
            array2d<real64, compflow::LAYOUT_PHASE> m_phaseMaxHistoricalVolFraction;

            //needed in hysteresis of wetting phase
            real64 m_phaseIntermediateMinVolFraction;

        };

        GEOS_HOST_DEVICE
        inline void TableCapillaryPressureHysteresis::KernelWrapper::compute(
                arraySlice1d<real64 const, compflow::USD_PHASE - 1> const &phaseVolFraction,
                arraySlice1d<real64 const, compflow::USD_PHASE - 1> const &phaseMaxHistoricalVolFraction,
                arraySlice1d<real64 const, compflow::USD_PHASE - 1> const &phaseMinHistoricalVolFraction,
                arraySlice1d<real64, cappres::USD_CAPPRES - 2> const &phaseTrappedVolFrac,
                arraySlice1d<real64, cappres::USD_CAPPRES - 2> const &phaseCapPressure,
                arraySlice2d<real64, cappres::USD_CAPPRES_DS - 2> const &dPhaseCapPressure_dPhaseVolFrac,
                ModeIndexType &mode
        ) const {
            LvArray::forValuesInSlice(dPhaseCapPressure_dPhaseVolFrac, [](real64 &val) { val = 0.0; });

            using PT = CapillaryPressureBase::PhaseType;
            integer const ipWater = m_phaseOrder[PT::WATER];
            integer const ipOil = m_phaseOrder[PT::OIL];
            integer const ipGas = m_phaseOrder[PT::GAS];

            if (ipWater >= 0 && ipOil >= 0 && ipGas >= 0) {
                computeThreePhase(ipWater, // wetting
                                  ipOil,   // intermediate
                                  ipGas,   // non-wetting
                                  phaseVolFraction,
                                  phaseMaxHistoricalVolFraction,
                                  phaseMinHistoricalVolFraction,
                                  phaseTrappedVolFrac,
                                  phaseCapPressure,
                                  dPhaseCapPressure_dPhaseVolFrac,
                                  mode);

            } else if (ipWater < 0) {
                computeTwoPhaseNonWetting(ipOil, // leading
                                          ipGas, // deduced
                                          phaseVolFraction,
                                          phaseMaxHistoricalVolFraction,
                                          phaseMinHistoricalVolFraction,
                                          phaseTrappedVolFrac,
                                          phaseCapPressure,
                                          dPhaseCapPressure_dPhaseVolFrac,
                                          mode);
            } else if (ipOil < 0) {
                computeTwoPhaseWetting(ipWater, // leading
                                       ipGas,   // deduced
                                       phaseVolFraction,
                                       phaseMaxHistoricalVolFraction,
                                       phaseMinHistoricalVolFraction,
                                       phaseTrappedVolFrac,
                                       phaseCapPressure,
                                       dPhaseCapPressure_dPhaseVolFrac,
                                       mode);
            } else if (ipGas < 0) {
                computeTwoPhaseWetting(ipWater, //leading
                                       ipOil,   //deduced
                                       phaseVolFraction,
                                       phaseMaxHistoricalVolFraction,
                                       phaseMinHistoricalVolFraction,
                                       phaseTrappedVolFrac,
                                       phaseCapPressure,
                                       dPhaseCapPressure_dPhaseVolFrac,
                                       mode);
            }


        }

        GEOS_HOST_DEVICE
        inline void TableCapillaryPressureHysteresis::KernelWrapper::update(const geos::localIndex k,
                                                                            const geos::localIndex q,
                                                                            const arraySlice1d<const geos::real64,
                                                                                    compflow::USD_PHASE
                                                                                    - 1> &phaseVolFraction) const {
            compute(phaseVolFraction,
                    m_phaseMaxHistoricalVolFraction[k],
                    m_phaseMinHistoricalVolFraction[k],
                    m_phaseTrappedVolFrac[k][q],
                    m_phaseCapPressure[k][q],
                    m_dPhaseCapPressure_dPhaseVolFrac[k][q],
                    m_mode[k]);
        }


    } //constitutive
} // geos

#endif //GEOS_CONSTITUTIVE_TABLECAPILLARYPRESSUREHYSTERESIS_HPP
