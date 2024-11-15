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

/**
 * @file MpiWrapper.hpp
 */

#ifndef GEOS_COMMON_MPIWRAPPER_HPP_
#define GEOS_COMMON_MPIWRAPPER_HPP_

#include "common/DataTypes.hpp"
#include "common/Span.hpp"

#if defined(GEOS_USE_MPI)
  #include <mpi.h>
#define MPI_PARAM( x ) x
#else
#define MPI_PARAM( x )
typedef int MPI_Comm;

#define MPI_COMM_NULL      ((MPI_Comm)0x04000000)
#define MPI_COMM_WORLD     ((MPI_Comm)0x44000000)
#define MPI_COMM_SELF      ((MPI_Comm)0x40000000)


typedef int MPI_Datatype;
#define MPI_CHAR           ((MPI_Datatype)0x4c000101)
#define MPI_SIGNED_CHAR    ((MPI_Datatype)0x4c000118)
#define MPI_UNSIGNED_CHAR  ((MPI_Datatype)0x4c000102)
#define MPI_BYTE           ((MPI_Datatype)0x4c00010d)
#define MPI_WCHAR          ((MPI_Datatype)0x4c00040e)
#define MPI_SHORT          ((MPI_Datatype)0x4c000203)
#define MPI_UNSIGNED_SHORT ((MPI_Datatype)0x4c000204)
#define MPI_INT            ((MPI_Datatype)0x4c000405)
#define MPI_UNSIGNED       ((MPI_Datatype)0x4c000406)
#define MPI_LONG           ((MPI_Datatype)0x4c000807)
#define MPI_UNSIGNED_LONG  ((MPI_Datatype)0x4c000808)
#define MPI_FLOAT          ((MPI_Datatype)0x4c00040a)
#define MPI_DOUBLE         ((MPI_Datatype)0x4c00080b)
#define MPI_LONG_DOUBLE    ((MPI_Datatype)0x4c00100c)
#define MPI_LONG_LONG_INT  ((MPI_Datatype)0x4c000809)
#define MPI_UNSIGNED_LONG_LONG ((MPI_Datatype)0x4c000819)
#define MPI_LONG_LONG      MPI_LONG_LONG_INT

typedef int MPI_Op;

#define MPI_MAX     (MPI_Op)(0x58000001)
#define MPI_MIN     (MPI_Op)(0x58000002)
#define MPI_SUM     (MPI_Op)(0x58000003)
#define MPI_PROD    (MPI_Op)(0x58000004)
#define MPI_LAND    (MPI_Op)(0x58000005)
#define MPI_BAND    (MPI_Op)(0x58000006)
#define MPI_LOR     (MPI_Op)(0x58000007)
#define MPI_BOR     (MPI_Op)(0x58000008)
#define MPI_LXOR    (MPI_Op)(0x58000009)
#define MPI_BXOR    (MPI_Op)(0x5800000a)
#define MPI_MINLOC  (MPI_Op)(0x5800000b)
#define MPI_MAXLOC  (MPI_Op)(0x5800000c)
#define MPI_REPLACE (MPI_Op)(0x5800000d)
#define MPI_NO_OP   (MPI_Op)(0x5800000e)

#define MPI_SUCCESS          0        /* Successful return code */
#define MPI_UNDEFINED (-32766)
#define MPI_STATUS_IGNORE (MPI_Status *)1
#define MPI_STATUSES_IGNORE (MPI_Status *)1
#define MPI_REQUEST_NULL  ((MPI_Request)0x2c000000)
typedef int MPI_Request;

typedef int MPI_Info;
#define MPI_INFO_NULL (MPI_Info)(0x60000000)

struct MPI_Status
{
  int junk;
};

#endif

#if defined(NDEBUG)
#define MPI_CHECK_ERROR( error ) ((void) error)
#else
#define MPI_CHECK_ERROR( error ) GEOS_ERROR_IF_NE( error, MPI_SUCCESS );
#endif


namespace geos
{

/// Global MPI communicator used by GEOSX.
#ifdef GEOS_USE_MPI
extern MPI_Comm MPI_COMM_GEOS;
#else
extern int MPI_COMM_GEOS;
#endif

/**
 * @struct MpiWrapper
 * This struct is a wrapper for all mpi.h functions that are used in GEOSX, and provides a collection of
 * convenience functions to make using the raw mpi functions simpler.
 *
 * The static wrapper functions around the mpi.h function are named by removing the "MPI_" from the beginning
 * of the native mpi function name. For instance the "Comm_rank()" function calls "MPI_Comm_rank()". Since all
 * wrapper functions are static, the should be referred to by their scoped name, for example "MpiWrapper::commRank()".
 */
struct MpiWrapper
{
public:

  /**
   * @enum Reduction
   * Strongly typed enum class for calling collective functions using MPI_Op
   */
  enum class Reduction
  {
    Max,  //!< Max
    Min,  //!< Min
    Sum,  //!< Sum
    Prod, //!< Prod
  };

  MpiWrapper() = delete;

  /**
   * @name FUNCTION GROUP for the direct wrappers around naitive MPI functions
   * @param[in] sendbuf  Pointer to the memory to read the sent data from.
   * @param[out] recvbuf Pointer to the memory to write the received data in.
   * @param[in] count    The number of data entries that are being communicated.
   * @param[in] datatype The MPI_Datatype that is being communicated.
   * @param[in] op       The collective MPI_Op to apply for the function.
   * @param[in] comm     The MPI_Comm communicator that the function will act on.
   *
   * Please see standard MPI documentation for a detailed description of the parameters for each function that
   * is being wrapped
   */
  ///@{

  static void barrier( MPI_Comm const & MPI_PARAM( comm )=MPI_COMM_GEOS );

  static int cartCoords( MPI_Comm comm, int rank, int maxdims, int coords[] );

  static int cartCreate( MPI_Comm comm_old, int ndims, const int dims[], const int periods[],
                         int reorder, MPI_Comm * comm_cart );

  static int cartRank( MPI_Comm comm, const int coords[] );

  static void commFree( MPI_Comm & comm );

  static int commRank( MPI_Comm const & MPI_PARAM( comm )=MPI_COMM_GEOS );

  static int commSize( MPI_Comm const & MPI_PARAM( comm )=MPI_COMM_GEOS );

  static bool commCompare( MPI_Comm const & comm1, MPI_Comm const & comm2 );

  static bool initialized();

  static int init( int * argc, char * * * argv );

  static void finalize();

  static MPI_Comm commDup( MPI_Comm const comm );

  static MPI_Comm commSplit( MPI_Comm const comm, int color, int key );

  static int test( MPI_Request * request, int * flag, MPI_Status * status );

  static int testAny( int count, MPI_Request array_of_requests[], int * idx, int * flags, MPI_Status array_of_statuses[] );

  static int testSome( int count, MPI_Request array_of_requests[], int * outcount, int array_of_indices[], MPI_Status array_of_statuses[] );

  static int testAll( int count, MPI_Request array_of_requests[], int * flags, MPI_Status array_of_statuses[] );

  /**
   * The same as test but doesn't deallocate requests regardless of their source.
   * @param[in] request The MPI_Request to check for completion
   * @param[out] flag Whether the request has completed or not
   * @param[out] status The current status of the request
   **/
  static int check( MPI_Request * request, int * flag, MPI_Status * status );

  /**
   * The same as testany but doesn't deallocate requests regardless of their source.
   * @note Since this doesn't deallocate or set request to MPI_REQUEST_NULL repeated calls
   *       to it with the same set of requests will return the first completed request each time
   *       thus if the first request is complete, that will always be returned unless the user
   *       deallocates/overwrites the request with MPI_REQUEST_NULL
   * @param[in] count The number of requests in the array to check
   * @param[in] request The MPI_Requests to check for completion
   * @param[out] idx The index of the first request in the array that is complete
   * @param[out] flag Whether a request has completed or not
   * @param[out] status The current status of all requests
   **/
  static int checkAny( int count, MPI_Request array_of_requests[], int * idx, int * flag, MPI_Status array_of_statuses[] );

  /**
   * The same as testall but doesn't deallocate requests regardless of their source.
   * @note Since this doesn't deallocate or set request to MPI_REQUEST_NULL repeated calls
   *       to it with a set of already-completed requests will continue to return.
   * @param[in] count The number of requests in the array to check
   * @param[in] request The MPI_Requests to check for completion
   * @param[out] flag Whether all requests have completed or not
   * @param[out] status The current status of all requests
   **/
  static int checkAll( int count, MPI_Request array_of_requests[], int * flag, MPI_Status array_of_statuses[] );

  static int wait( MPI_Request * request, MPI_Status * status );

  static int waitAny( int count, MPI_Request array_of_requests[], int * indx, MPI_Status array_of_statuses[] );

  static int waitSome( int count, MPI_Request array_of_requests[], int * outcount, int array_of_indices[], MPI_Status array_of_statuses[] );

  static int waitAll( int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[] );

  static double wtime( void );


  /**
   * Wait on MPI_Requests to complete on at a time and trigger a callback to
   *  process the completion.
   * @param[in] count The number of MPI_Requests being processed.
   * @param[inout] array_of_requests The MPI_Requests to actively wait on.
   * @param[in] func A callable object accepting an integer denoting the MPI_Request index
   *              which has completed.
   * @return MPI_SUCCESS or an MPI_ERROR returned by internal calls to MPI_WaitAny.
   */
  static int activeWaitAny( const int count,
                            MPI_Request array_of_requests[],
                            MPI_Status array_of_statuses[],
                            std::function< MPI_Request ( int ) > func );

  /**
   * Wait on MPI_Requests to complete on or more at a time and trigger a callback to
   *  process the completion.
   * @param[in] count The number of MPI_Requests being processed.
   * @param[inout] array_of_requests The MPI_Requests to actively wait on.
   * @param[in] func A callable object accepting an integer denoting the MPI_Request index
   *              which has completed.
   * @return MPI_SUCCESS or an MPI_ERROR returned by internal calls to MPI_WaitSome.
   */
  static int activeWaitSome( const int count,
                             MPI_Request array_of_requests[],
                             MPI_Status array_of_statuses[],
                             std::function< MPI_Request ( int ) > func );

  /**
   * Active non-blocking phased communication with multiple participants,
   *  each participant in each phase may depend on the previous phases
   *  being fully complete prior to entry into a subsequent phase.
   * @param[in] participants The number of participants in each phase
   * @param[in] phases A vector of function objects taking int and returning MPI_Requests
   *               denoting the state of that participant in that phase.
   * @note The restriction on phase[N](index) being called is that phase[N-1](0 - (particpants-1))
   *       have all been called and the MPI_Requests from those calls have all been completed.
   * @note One can add a final recv phase by having that phase return MPI_REQUEST_NULL.
   * @return MPI_SUCCESS or and MPI_ERROR from internal calls to MPI_WaitAny.
   */
  static int activeWaitSomeCompletePhase( const int participants,
                                          std::vector< std::tuple< MPI_Request *, MPI_Status *, std::function< MPI_Request ( int ) > > > const & phases );

  /**
   * Active blocking phased communication with multiple participants,
   *  each participant in each phase may depend on the previous phases
   *  being fully complete prior to entry into a subsequent phase.
   * @param[in] participants The number of participants in each phase
   * @param[in] phases A vector of function objects taking int and returning MPI_Requests
   *               denoting the state of that participant in that phase.
   * @note The restriction on phase[N](index) being called is that phase[N-1](0 - (particpants-1))
   *       and phase[N](0..index) have all been called and the MPI_Requests from those calls have all
   *       been completed.
   * @note One can add a final recv phase by having that phase return MPI_REQUEST_NULL.
   * @return MPI_SUCCESS or and MPI_ERROR from internal calls to MPI_WaitAny.
   */
  static int activeWaitOrderedCompletePhase( const int participants,
                                             std::vector< std::tuple< MPI_Request *, MPI_Status *, std::function< MPI_Request ( int ) > > > const & phases );
  ///@}

#if !defined(GEOS_USE_MPI)
  static std::map< int, std::pair< int, void * > > & getTagToPointersMap()
  {
    static std::map< int, std::pair< int, void * > > tagToPointers;
    return tagToPointers;
  }
#endif

  /**
   * @brief Compute the number of ranks allocated on the same node
   * @return The number of MPI ranks on the current node.
   */
  static int nodeCommSize();

  /**
   * @brief Strongly typed wrapper around MPI_Allgather.
   * @tparam T_SEND The pointer type for \p sendbuf
   * @tparam T_RECV The pointer type for \p recvbuf
   * @param[in] sendbuf The pointer to the sending buffer.
   * @param[in] sendcount The number of values to send.
   * @param[out] recvbuf The pointer to the receive buffer.
   * @param[in] recvcount The number of values to receive.
   * @param[in] comm The MPI_Comm over which the gather operates.
   * @return The return value of the underlying call to MPI_Allgather().
   */
  template< typename T_SEND, typename T_RECV >
  static int allgather( T_SEND const * sendbuf,
                        int sendcount,
                        T_RECV * recvbuf,
                        int recvcount,
                        MPI_Comm comm );

  /**
   * @brief Strongly typed wrapper around MPI_Allgatherv.
   * @tparam T_SEND The pointer type for \p sendbuf
   * @tparam T_RECV The pointer type for \p recvbuf
   * @param[in] sendbuf The pointer to the sending buffer.
   * @param[in] sendcount The number of values to send.
   * @param[out] recvbuf The pointer to the receive buffer.
   * @param[in] recvcounts The number of values to receive.
   * @param[in] displacements An array containing the displacement to apply to the message received by each process
   * @param[in] comm The MPI_Comm over which the gather operates.
   * @return The return value of the underlying call to MPI_Allgatherv().
   */
  template< typename T_SEND, typename T_RECV >
  static int allgatherv( T_SEND const * sendbuf,
                         int sendcount,
                         T_RECV * recvbuf,
                         int * recvcounts,
                         int * displacements,
                         MPI_Comm comm );

  /**
   * @brief Convenience function for MPI_Allgather.
   * @tparam T The type to send/recieve. This must have a valid conversion to MPI_Datatype in getMpiType();
   * @param[in] myValue The value to send.
   * @param[out] allValues The values recived from each rank.
   */
  template< typename T >
  static void allGather( T const myValue, array1d< T > & allValues, MPI_Comm comm = MPI_COMM_GEOS );

  template< typename T >
  static int allGather( arrayView1d< T const > const & sendbuf,
                        array1d< T > & recvbuf,
                        MPI_Comm comm = MPI_COMM_GEOS );

  /**
   * @brief Strongly typed wrapper around MPI_Allreduce.
   * @param[in] sendbuf The pointer to the sending buffer.
   * @param[out] recvbuf The pointer to the receive buffer.
   * @param[in] count The number of values to send/receive.
   * @param[in] op The MPI_Op to perform.
   * @param[in] comm The MPI_Comm over which the gather operates.
   * @return The return value of the underlying call to MPI_Allreduce().
   */
  template< typename T >
  static int allReduce( T const * sendbuf, T * recvbuf, int count, MPI_Op op, MPI_Comm comm = MPI_COMM_GEOS );

  /**
   * @brief Convenience wrapper for the MPI_Allreduce function.
   * @tparam T type of data to reduce. Must correspond to a valid MPI_Datatype.
   * @param value The value to send to the reduction.
   * @param op The Reduction enum to perform.
   * @param comm The communicator.
   * @return The value of reduction across all ranks
   */
  template< typename T >
  static T allReduce( T const & value, Reduction const op, MPI_Comm comm = MPI_COMM_GEOS );

  /**
   * @brief Convenience wrapper for the MPI_Allreduce function. Version for sequences.
   * @tparam T type of data to reduce. Must correspond to a valid MPI_Datatype.
   * @param src[in] The values to send to the reduction.
   * @param dst[out] The resulting values.
   * @param op The Reduction enum to perform.
   * @param comm The communicator.
   */
  template< typename T >
  static void allReduce( Span< T const > src, Span< T > dst, Reduction const op, MPI_Comm comm = MPI_COMM_GEOS );


  /**
   * @brief Strongly typed wrapper around MPI_Reduce.
   * @param[in] sendbuf The pointer to the sending buffer.
   * @param[out] recvbuf The pointer to the receive buffer (only significant at root).
   * @param[in] count The number of values to send/receive.
   * @param[in] op The MPI_Op to perform.
   * @param[in] comm The MPI_Comm over which the gather operates.
   * @return The return value of the underlying call to MPI_Reduce().
   */
  template< typename T >
  static int reduce( T const * sendbuf, T * recvbuf, int count, MPI_Op op, int root, MPI_Comm comm = MPI_COMM_GEOS );

  /**
   * @brief Convenience wrapper for the MPI_Reduce function.
   * @tparam T type of data to reduce. Must correspond to a valid MPI_Datatype.
   * @param value The value to send to the reduction.
   * @param op The Reduction enum to perform.
   * @param comm The communicator.
   * @return The value of reduction (only significant at root)
   */
  template< typename T >
  static T reduce( T const & value, Reduction const op, int root, MPI_Comm comm = MPI_COMM_GEOS );

  /**
   * @brief Convenience wrapper for the MPI_Reduce function. Version for sequences.
   * @tparam T type of data to reduce. Must correspond to a valid MPI_Datatype.
   * @param src[in] The values to send to the reduction.
   * @param dst[out] The resulting values (only significant at root).
   * @param op The Reduction enum to perform.
   * @param comm The communicator.
   */
  template< typename T >
  static void reduce( Span< T const > src, Span< T > dst, Reduction const op, int root, MPI_Comm comm = MPI_COMM_GEOS );


  template< typename T >
  static int scan( T const * sendbuf, T * recvbuf, int count, MPI_Op op, MPI_Comm comm );

  template< typename T >
  static int exscan( T const * sendbuf, T * recvbuf, int count, MPI_Op op, MPI_Comm comm );

  /**
   * @brief Strongly typed wrapper around MPI_Bcast.
   * @param[in/out] buffer The pointer to the send/recv buffer
   * @param[in] count The number of data types to send.
   * @param[in] root The rank sending the data.
   * @param[in] comm The MPI_Comm over which the MPI_Bcast operates.
   * @return The return value of the underlying call to MPI_Bcast().
   */
  template< typename T >
  static int bcast( T * buffer, int count, int root, MPI_Comm comm );


  /**
   * @brief Convenience function for MPI_Broadcast.
   * @tparam T The type to send/recieve. This must have a valid conversion to MPI_Type in getMpiType();
   * @param[in/out] myValue The value to send from the \p srcRank to the other ranks.
   * @param srcRank The rank that is sending the \p value.
   */
  template< typename T >
  static void broadcast( T & value, int srcRank = 0, MPI_Comm comm = MPI_COMM_GEOS );

  /**
   * @brief Strongly typed wrapper around MPI_Gather().
   * @tparam TS The pointer type for \p sendbuf
   * @tparam TR The pointer type for \p recvbuf
   * @param[in] sendbuf The pointer to the sending buffer.
   * @param[in] sendcount The number of values to send.
   * @param[out] recvbuf The pointer to the receive buffer.
   * @param[in] recvcount The number of values to receive.
   * @param[in] root The rank recieving the data.
   * @param[in] comm The MPI_Comm over which the gather operates.
   * @return
   */
  template< typename TS, typename TR >
  static int gather( TS const * const sendbuf,
                     int sendcount,
                     TR * const recvbuf,
                     int recvcount,
                     int root,
                     MPI_Comm comm );

  /**
   * @brief Strongly typed wrapper around MPI_Gatherv.
   * @tparam TS The pointer type for \p sendbuf
   * @tparam TR The pointer type for \p recvbuf
   * @param[in] sendbuf The pointer to the sending buffer.
   * @param[in] sendcount The number of values to send.
   * @param[out] recvbuf The pointer to the receive buffer.
   * @param[in] recvcount The number of values to receive.
   * @param[in] displs integer array (of length group size). Entry i specifies the displacement relative to recvbuf at
   *                   which to place the incoming data from process i (significant only at root).
   * @param[in] root The rank recieving the data.
   * @param[in] comm The MPI_Comm over which the gather operates.
   * @return
   */
  template< typename TS, typename TR >
  static int gatherv( TS const * const sendbuf,
                      int sendcount,
                      TR * const recvbuf,
                      const int * recvcounts,
                      const int * displs,
                      int root,
                      MPI_Comm comm );

  /**
   * @brief Returns an MPI_Op associated with our strongly typed Reduction enum.
   * @param[in] op The value of the Reduction enum to get an MPI_Op for.
   * @return The MPI_Op associated with \p op.
   */
  static MPI_Op getMpiOp( Reduction const op );

  template< typename T >
  static int recv( array1d< T > & buf,
                   int MPI_PARAM( source ),
                   int tag,
                   MPI_Comm MPI_PARAM( comm ),
                   MPI_Status * MPI_PARAM( request ) );

  template< typename T >
  static int iSend( arrayView1d< T > const & buf,
                    int MPI_PARAM( dest ),
                    int tag,
                    MPI_Comm MPI_PARAM( comm ),
                    MPI_Request * MPI_PARAM( request ) );

  /**
   * @brief Strongly typed wrapper around MPI_Irecv()
   * @param[out] buf The pointer to the buffer that contains the data to be received.
   * @param[in] count The number of elements in \p buf
   * @param[in] source The rank of the source process within \p comm.
   * @param[in] tag The message tag that is be used to distinguish different types of messages
   * @param[in] comm The handle to the MPI_Comm
   * @param[out] request Pointer to the MPI_Request associated with this request.
   * @return
   */
  template< typename T >
  static int iRecv( T * const buf,
                    int count,
                    int source,
                    int tag,
                    MPI_Comm comm,
                    MPI_Request * request );

  /**
   * @brief Strongly typed wrapper around MPI_Isend()
   * @param[in] buf The pointer to the buffer that contains the data to be sent.
   * @param[in] count The number of elements in \p buf.
   * @param[in] dest The rank of the destination process within \p comm.
   * @param[in] tag The message tag that is be used to distinguish different types of messages.
   * @param[in] comm The handle to the MPI_Comm.
   * @param[out] request Pointer to the MPI_Request associated with this request.
   * @return
   */
  template< typename T >
  static int iSend( T const * const buf,
                    int count,
                    int dest,
                    int tag,
                    MPI_Comm comm,
                    MPI_Request * request );

  /**
   * @brief Compute exclusive prefix sum and full sum
   * @tparam T type of local (rank) value
   * @tparam U type of global (sum) value
   * @param[in] value the local value
   * @return a pair where first is the prefix sum, second is the full sum
   */
  template< typename U, typename T >
  static U prefixSum( T const value, MPI_Comm comm = MPI_COMM_GEOS );

  /**
   * @brief Convenience function for a MPI_Allreduce using a MPI_SUM operation.
   * @param[in] value the value to send into the reduction.
   * @return The sum of all \p value across the ranks.
   */
  template< typename T >
  static T sum( T const & value, MPI_Comm comm = MPI_COMM_GEOS );

  /**
   * @brief Convenience function for a MPI_Allreduce using a MPI_SUM operation.
   * @param[in] src the value to send into the reduction.
   * @param[out] dst The resulting values.
   * @return The sum of all \p value across the ranks.
   */
  template< typename T >
  static void sum( Span< T const > src, Span< T > dst, MPI_Comm comm = MPI_COMM_GEOS );

  /**
   * @brief Convenience function for a MPI_Allreduce using a MPI_MIN operation.
   * @param value the value to send into the reduction.
   * @return The minimum of all \p value across the ranks.
   */
  template< typename T >
  static T min( T const & value, MPI_Comm comm = MPI_COMM_GEOS );

  /**
   * @brief Convenience function for a MPI_Allreduce using a MPI_MIN operation.
   * @param[in] src the value to send into the reduction.
   * @param[out] dst The resulting values.
   * @return The minimum of all \p value across the ranks.
   */
  template< typename T >
  static void min( Span< T const > src, Span< T > dst, MPI_Comm comm = MPI_COMM_GEOS );

  /**
   * @brief Convenience function for a MPI_Allreduce using a MPI_MAX operation.
   * @param[in] value the value to send into the reduction.
   * @return The maximum of all \p value across the ranks.
   */
  template< typename T >
  static T max( T const & value, MPI_Comm comm = MPI_COMM_GEOS );

  /**
   * @brief Convenience function for a MPI_Allreduce using a MPI_MAX operation.
   * @param[in] value the value to send into the reduction.
   * @param[out] dst The resulting values.
   * @return The maximum of all \p value across the ranks.
   */
  template< typename T >
  static void max( Span< T const > src, Span< T > dst, MPI_Comm comm = MPI_COMM_GEOS );


  /**
   * @brief Convenience function for MPI_Gather using a MPI_MAX operation on struct of value and location
   * @brief Max is performed on value and location (global index) is returned
   * @param[in] struct to send into the max gather.
   * @return struct with max val and location
   */
  template< typename T > static T maxValLoc( T localValueLocation, MPI_Comm comm = MPI_COMM_GEOS );

};

namespace internal
{

template< typename T, typename ENABLE = void >
struct MpiTypeImpl {};

#define ADD_MPI_TYPE_MAP( T, MPI_T ) \
  template<> struct MpiTypeImpl< T > { static MPI_Datatype get() { return MPI_T; } }

ADD_MPI_TYPE_MAP( float, MPI_FLOAT );
ADD_MPI_TYPE_MAP( double, MPI_DOUBLE );

ADD_MPI_TYPE_MAP( char, MPI_CHAR );
ADD_MPI_TYPE_MAP( signed char, MPI_SIGNED_CHAR );
ADD_MPI_TYPE_MAP( unsigned char, MPI_UNSIGNED_CHAR );

ADD_MPI_TYPE_MAP( int, MPI_INT );
ADD_MPI_TYPE_MAP( long int, MPI_LONG );
ADD_MPI_TYPE_MAP( long long int, MPI_LONG_LONG );

ADD_MPI_TYPE_MAP( unsigned int, MPI_UNSIGNED );
ADD_MPI_TYPE_MAP( unsigned long int, MPI_UNSIGNED_LONG );
ADD_MPI_TYPE_MAP( unsigned long long int, MPI_UNSIGNED_LONG_LONG );

#undef ADD_MPI_TYPE_MAP

template< typename T >
struct MpiTypeImpl< T, std::enable_if_t< std::is_enum< T >::value > >
{
  static MPI_Datatype get() { return MpiTypeImpl< std::underlying_type_t< T > >::get(); }
};

template< typename T >
MPI_Datatype getMpiType()
{
  return MpiTypeImpl< T >::get();
}

}

inline MPI_Op MpiWrapper::getMpiOp( Reduction const op )
{
  switch( op )
  {
    case Reduction::Sum:
    {
      return MPI_SUM;
    }
    case Reduction::Min:
    {
      return MPI_MIN;
    }
    case Reduction::Max:
    {
      return MPI_MAX;
    }
    case Reduction::Prod:
    {
      return MPI_PROD;
    }
    default:
      GEOS_ERROR( "Unsupported reduction operation" );
      return MPI_NO_OP;
  }
}

template< typename T_SEND, typename T_RECV >
int MpiWrapper::allgather( T_SEND const * const sendbuf,
                           int sendcount,
                           T_RECV * const recvbuf,
                           int recvcount,
                           MPI_Comm MPI_PARAM( comm ) )
{
#ifdef GEOS_USE_MPI
  return MPI_Allgather( sendbuf, sendcount, internal::getMpiType< T_SEND >(),
                        recvbuf, recvcount, internal::getMpiType< T_RECV >(),
                        comm );
#else
  static_assert( std::is_same< T_SEND, T_RECV >::value,
                 "MpiWrapper::allgather() for serial run requires send and receive buffers are of the same type" );
  GEOS_ERROR_IF_NE_MSG( sendcount, recvcount, "sendcount is not equal to recvcount." );
  std::copy( sendbuf, sendbuf + sendcount, recvbuf )
  return 0;
#endif
}

template< typename T_SEND, typename T_RECV >
int MpiWrapper::allgatherv( T_SEND const * const sendbuf,
                            int sendcount,
                            T_RECV * const recvbuf,
                            int * recvcounts,
                            int * displacements,
                            MPI_Comm MPI_PARAM( comm ) )
{
#ifdef GEOS_USE_MPI
  return MPI_Allgatherv( sendbuf, sendcount, internal::getMpiType< T_SEND >(),
                         recvbuf, recvcounts, displacements, internal::getMpiType< T_RECV >(),
                         comm );
#else
  static_assert( std::is_same< T_SEND, T_RECV >::value,
                 "MpiWrapper::allgatherv() for serial run requires send and receive buffers are of the same type" );
  GEOS_ERROR_IF_NE_MSG( sendcount, recvcount, "sendcount is not equal to recvcount." );
  std::copy( sendbuf, sendbuf + sendcount, recvbuf )
  return 0;
#endif
}


template< typename T >
void MpiWrapper::allGather( T const myValue, array1d< T > & allValues, MPI_Comm MPI_PARAM( comm ) )
{
#ifdef GEOS_USE_MPI
  int const mpiSize = commSize( comm );
  allValues.resize( mpiSize );

  MPI_Datatype const MPI_TYPE = internal::getMpiType< T >();

  MPI_Allgather( &myValue, 1, MPI_TYPE, allValues.data(), 1, MPI_TYPE, comm );

#else
  allValues.resize( 1 );
  allValues[0] = myValue;
#endif
}

template< typename T >
int MpiWrapper::allGather( arrayView1d< T const > const & sendValues,
                           array1d< T > & allValues,
                           MPI_Comm MPI_PARAM( comm ) )
{
  int const sendSize = LvArray::integerConversion< int >( sendValues.size() );
#ifdef GEOS_USE_MPI
  int const mpiSize = commSize( comm );
  allValues.resize( mpiSize * sendSize );
  return MPI_Allgather( sendValues.data(),
                        sendSize,
                        internal::getMpiType< T >(),
                        allValues.data(),
                        sendSize,
                        internal::getMpiType< T >(),
                        comm );

#else
  allValues.resize( sendSize );
  for( localIndex a=0; a<sendSize; ++a )
  {
    allValues[a] = sendValues[a];
  }
  return 0;
#endif
}

template< typename T >
int MpiWrapper::allReduce( T const * const sendbuf,
                           T * const recvbuf,
                           int const count,
                           MPI_Op const MPI_PARAM( op ),
                           MPI_Comm const MPI_PARAM( comm ) )
{
#ifdef GEOS_USE_MPI
  MPI_Datatype const mpiType = internal::getMpiType< T >();
  return MPI_Allreduce( sendbuf == recvbuf ? MPI_IN_PLACE : sendbuf, recvbuf, count, mpiType, op, comm );
#else
  if( sendbuf != recvbuf )
  {
    memcpy( recvbuf, sendbuf, count * sizeof( T ) );
  }
  return 0;
#endif
}

template< typename T >
int MpiWrapper::reduce( T const * const sendbuf,
                        T * const recvbuf,
                        int const count,
                        MPI_Op const MPI_PARAM( op ),
                        int root,
                        MPI_Comm const MPI_PARAM( comm ) )
{
#ifdef GEOS_USE_MPI
  MPI_Datatype const mpiType = internal::getMpiType< T >();
  return MPI_Reduce( sendbuf == recvbuf ? MPI_IN_PLACE : sendbuf, recvbuf, count, mpiType, op, root, comm );
#else
  if( sendbuf != recvbuf )
  {
    memcpy( recvbuf, sendbuf, count * sizeof( T ) );
  }
  return 0;
#endif
}

template< typename T >
int MpiWrapper::scan( T const * const sendbuf,
                      T * const recvbuf,
                      int count,
                      MPI_Op MPI_PARAM( op ),
                      MPI_Comm MPI_PARAM( comm ) )
{
#ifdef GEOS_USE_MPI
  return MPI_Scan( sendbuf, recvbuf, count, internal::getMpiType< T >(), op, comm );
#else
  memcpy( recvbuf, sendbuf, count*sizeof(T) );
  return 0;
#endif
}

template< typename T >
int MpiWrapper::exscan( T const * const MPI_PARAM( sendbuf ),
                        T * const recvbuf,
                        int count,
                        MPI_Op MPI_PARAM( op ),
                        MPI_Comm MPI_PARAM( comm ) )
{
#ifdef GEOS_USE_MPI
  return MPI_Exscan( sendbuf, recvbuf, count, internal::getMpiType< T >(), op, comm );
#else
  memset( recvbuf, 0, count*sizeof(T) );
  return 0;
#endif
}

template< typename T >
int MpiWrapper::bcast( T * const MPI_PARAM( buffer ),
                       int MPI_PARAM( count ),
                       int MPI_PARAM( root ),
                       MPI_Comm MPI_PARAM( comm ) )
{
#ifdef GEOS_USE_MPI
  return MPI_Bcast( buffer, count, internal::getMpiType< T >(), root, comm );
#else
  return 0;
#endif

}

template< typename T >
void MpiWrapper::broadcast( T & MPI_PARAM( value ), int MPI_PARAM( srcRank ), MPI_Comm MPI_PARAM( comm ) )
{
#ifdef GEOS_USE_MPI
  MPI_Bcast( &value, 1, internal::getMpiType< T >(), srcRank, comm );
#endif
}

template<>
inline
void MpiWrapper::broadcast< string >( string & MPI_PARAM( value ),
                                      int MPI_PARAM( srcRank ),
                                      MPI_Comm MPI_PARAM( comm ) )
{
#ifdef GEOS_USE_MPI
  int size = LvArray::integerConversion< int >( value.size() );
  broadcast( size, srcRank, comm );
  value.resize( size );
  MPI_Bcast( const_cast< char * >( value.data() ), size, internal::getMpiType< char >(), srcRank, comm );
#endif
}

template< typename TS, typename TR >
int MpiWrapper::gather( TS const * const sendbuf,
                        int sendcount,
                        TR * const recvbuf,
                        int recvcount,
                        int MPI_PARAM( root ),
                        MPI_Comm MPI_PARAM( comm ) )
{
#ifdef GEOS_USE_MPI
  return MPI_Gather( sendbuf, sendcount, internal::getMpiType< TS >(),
                     recvbuf, recvcount, internal::getMpiType< TR >(),
                     root, comm );
#else
  static_assert( std::is_same< TS, TR >::value,
                 "MpiWrapper::gather() for serial run requires send and receive buffers are of the same type" );
  std::size_t const sendBufferSize = sendcount * sizeof(TS);
  std::size_t const recvBufferSize = recvcount * sizeof(TR);
  GEOS_ERROR_IF_NE_MSG( sendBufferSize, recvBufferSize, "size of send buffer and receive buffer are not equal" );
  memcpy( recvbuf, sendbuf, sendBufferSize );
  return 0;
#endif
}

template< typename TS, typename TR >
int MpiWrapper::gatherv( TS const * const sendbuf,
                         int sendcount,
                         TR * const recvbuf,
                         const int * recvcounts,
                         const int * MPI_PARAM( displs ),
                         int MPI_PARAM( root ),
                         MPI_Comm MPI_PARAM( comm ) )
{
#ifdef GEOS_USE_MPI
  return MPI_Gatherv( sendbuf, sendcount, internal::getMpiType< TS >(),
                      recvbuf, recvcounts, displs, internal::getMpiType< TR >(),
                      root, comm );
#else
  static_assert( std::is_same< TS, TR >::value,
                 "MpiWrapper::gather() for serial run requires send and receive buffers are of the same type" );
  std::size_t const sendBufferSize = sendcount * sizeof(TS);
  std::size_t const recvBufferSize = recvcounts[0] * sizeof(TR);
  GEOS_ERROR_IF_NE_MSG( sendBufferSize, recvBufferSize, "size of send buffer and receive buffer are not equal" );
  memcpy( recvbuf, sendbuf, sendBufferSize );
  return 0;
#endif
}

template< typename T >
int MpiWrapper::iRecv( T * const buf,
                       int count,
                       int MPI_PARAM( source ),
                       int tag,
                       MPI_Comm MPI_PARAM( comm ),
                       MPI_Request * MPI_PARAM( request ) )
{
#ifdef GEOS_USE_MPI
  GEOS_ERROR_IF( (*request)!=MPI_REQUEST_NULL,
                 "Attempting to use an MPI_Request that is still in use." );
  return MPI_Irecv( buf, count, internal::getMpiType< T >(), source, tag, comm, request );
#else
  std::map< int, std::pair< int, void * > > & pointerMap = getTagToPointersMap();
  std::map< int, std::pair< int, void * > >::iterator iPointer = pointerMap.find( tag );

  if( iPointer==pointerMap.end() )
  {
    pointerMap.insert( {tag, {1, buf} } );
  }
  else
  {
    GEOS_ERROR_IF( iPointer->second.first != 0,
                   "Tag does is assigned, but pointer was not set by iSend." );
    memcpy( buf, iPointer->second.second, count*sizeof(T) );
    pointerMap.erase( iPointer );
  }
  return 0;
#endif
}

template< typename T >
int MpiWrapper::recv( array1d< T > & buf,
                      int MPI_PARAM( source ),
                      int tag,
                      MPI_Comm MPI_PARAM( comm ),
                      MPI_Status * MPI_PARAM( request ) )
{
#ifdef GEOS_USE_MPI
  MPI_Status status;
  int count;
  MPI_Probe( source, tag, comm, &status );
  MPI_Get_count( &status, MPI_CHAR, &count );

  GEOS_ASSERT_EQ( count % sizeof( T ), 0 );
  buf.resize( count / sizeof( T ) );

  return MPI_Recv( reinterpret_cast< char * >( buf.data() ),
                   count,
                   MPI_CHAR,
                   source,
                   tag,
                   comm,
                   request );
#else
  GEOS_ERROR( "Not implemented!" );
  return MPI_SUCCESS;
#endif
}

template< typename T >
int MpiWrapper::iSend( arrayView1d< T > const & buf,
                       int MPI_PARAM( dest ),
                       int tag,
                       MPI_Comm MPI_PARAM( comm ),
                       MPI_Request * MPI_PARAM( request ) )
{
#ifdef GEOS_USE_MPI
  GEOS_ERROR_IF( (*request)!=MPI_REQUEST_NULL,
                 "Attempting to use an MPI_Request that is still in use." );
  return MPI_Isend( reinterpret_cast< void const * >( buf.data() ),
                    buf.size() * sizeof( T ),
                    MPI_CHAR,
                    dest,
                    tag,
                    comm,
                    request );
#else
  GEOS_ERROR( "Not implemented." );
  return MPI_SUCCESS;
#endif
}

template< typename T >
int MpiWrapper::iSend( T const * const buf,
                       int count,
                       int MPI_PARAM( dest ),
                       int tag,
                       MPI_Comm MPI_PARAM( comm ),
                       MPI_Request * MPI_PARAM( request ) )
{
#ifdef GEOS_USE_MPI
  GEOS_ERROR_IF( (*request)!=MPI_REQUEST_NULL,
                 "Attempting to use an MPI_Request that is still in use." );
  return MPI_Isend( buf, count, internal::getMpiType< T >(), dest, tag, comm, request );
#else
  std::map< int, std::pair< int, void * > > & pointerMap = getTagToPointersMap();
  std::map< int, std::pair< int, void * > >::iterator iPointer = pointerMap.find( tag );

  if( iPointer==pointerMap.end() )
  {
    pointerMap.insert( {tag, {0, const_cast< T * >(buf)}
                       } );
  }
  else
  {
    GEOS_ERROR_IF( iPointer->second.first != 1,
                   "Tag does is assigned, but pointer was not set by iRecv." );
    memcpy( iPointer->second.second, buf, count*sizeof(T) );
    pointerMap.erase( iPointer );
  }
  return 0;
#endif
}

template< typename U, typename T >
U MpiWrapper::prefixSum( T const value, MPI_Comm comm )
{
  U localResult;

#ifdef GEOS_USE_MPI
  U const convertedValue = value;
  int const error = MPI_Exscan( &convertedValue, &localResult, 1, internal::getMpiType< U >(), MPI_SUM, comm );
  MPI_CHECK_ERROR( error );
#endif
  if( commRank() == 0 )
  {
    localResult = 0;
  }

  return localResult;
}


template< typename T >
T MpiWrapper::allReduce( T const & value, Reduction const op, MPI_Comm const comm )
{
  T result;
  allReduce( &value, &result, 1, getMpiOp( op ), comm );
  return result;
}

template< typename T >
void MpiWrapper::allReduce( Span< T const > const src, Span< T > const dst, Reduction const op, MPI_Comm const comm )
{
  GEOS_ASSERT_EQ( src.size(), dst.size() );
  allReduce( src.data(), dst.data(), LvArray::integerConversion< int >( src.size() ), getMpiOp( op ), comm );
}

template< typename T >
T MpiWrapper::sum( T const & value, MPI_Comm comm )
{
  return MpiWrapper::allReduce( value, Reduction::Sum, comm );
}

template< typename T >
void MpiWrapper::sum( Span< T const > src, Span< T > dst, MPI_Comm comm )
{
  MpiWrapper::allReduce( src, dst, Reduction::Sum, comm );
}

template< typename T >
T MpiWrapper::min( T const & value, MPI_Comm comm )
{
  return MpiWrapper::allReduce( value, Reduction::Min, comm );
}

template< typename T >
void MpiWrapper::min( Span< T const > src, Span< T > dst, MPI_Comm comm )
{
  MpiWrapper::allReduce( src, dst, Reduction::Min, comm );
}

template< typename T >
T MpiWrapper::max( T const & value, MPI_Comm comm )
{
  return MpiWrapper::allReduce( value, Reduction::Max, comm );
}

template< typename T >
void MpiWrapper::max( Span< T const > src, Span< T > dst, MPI_Comm comm )
{
  MpiWrapper::allReduce( src, dst, Reduction::Max, comm );
}


template< typename T >
T MpiWrapper::reduce( T const & value, Reduction const op, int root, MPI_Comm const comm )
{
  T result;
  reduce( &value, &result, 1, getMpiOp( op ), root, comm );
  return result;
}

template< typename T >
void MpiWrapper::reduce( Span< T const > const src, Span< T > const dst, Reduction const op, int root, MPI_Comm const comm )
{
  GEOS_ASSERT_EQ( src.size(), dst.size() );
  reduce( src.data(), dst.data(), LvArray::integerConversion< int >( src.size() ), getMpiOp( op ), root, comm );
}

// Mpi helper function to return  struct containing the max value and location across ranks
template< typename T >
T MpiWrapper::maxValLoc( T localValueLocation, MPI_Comm comm )
{
  // Ensure T is trivially copyable
  static_assert( std::is_trivially_copyable< T >::value, "maxValLoc requires a trivially copyable type" );

  // T to have only 2 data members named value and location
  static_assert( (sizeof(T::value)+sizeof(T::location)) == sizeof(T) );

  // Ensure T has value and location members are scalars
  static_assert( std::is_scalar_v< decltype(T::value) > || std::is_scalar_v< decltype(T::location) >, "members of struct should be scalar" );
  static_assert( !std::is_pointer_v< decltype(T::value) > && !std::is_pointer_v< decltype(T::location) >, "members of struct should not be pointers" );

  // receive "buffer"
  int const numProcs =  commSize( comm );
  std::vector< T > recvValLoc( numProcs );

  MPI_Allgather( &localValueLocation, sizeof(T), MPI_BYTE, recvValLoc.data(), sizeof(T), MPI_BYTE, comm );

  T maxValLoc= *std::max_element( recvValLoc.begin(),
                                  recvValLoc.end(),
                                  []( auto & lhs, auto & rhs ) -> bool {return lhs.value  <  rhs.value; } );

  return maxValLoc;
}
} /* namespace geos */

#endif /* GEOS_COMMON_MPIWRAPPER_HPP_ */
