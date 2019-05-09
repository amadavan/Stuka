if (GUROBI_FOUND)
    set(GUROBI_SOURCES
            ../include/stuka/LP/gurobi.h
            ../include/stuka/LP/gurobi_lp.h
            ../include/stuka/QP/gurobi.h
            ../include/stuka/QP/gurobi_qp.h
            LP/gurobi.cpp
            LP/gurobi_lp.cpp
            QP/gurobi.cpp
            QP/gurobi_qp.cpp)
else ()
    set(GUROBI_LIBRARIES "")
endif ()

if (HDF5_FOUND)
    set(HDF5_SOURCES
            ../include/stuka/util/callback/save_hdf5.h
            util/callback/save_hdf5.cpp)
else ()
    set(HDF5_LIBRARIES "")
endif ()

set(STUKA_SOURCES
        ../include/stuka/LP/base_lp.h
        ../include/stuka/LP/base_solver.h
        ../include/stuka/LP/lp.h
        ../include/stuka/LP/mehrotra_pc.h
        ../include/stuka/LP/slack_lp.h
        ../include/stuka/QP/base_qp.h
        ../include/stuka/QP/base_solver.h
        ../include/stuka/QP/qp.h
        ../include/stuka/dLP/decomposed_lp.h
        ../include/stuka/dLP/base_solver.h
        ../include/stuka/dLP/subproblem.h
        ../include/stuka/dLP/cre.h
        ../include/stuka/dLP/cre_subproblem.h
        ../include/stuka/dLP/critical_region.h
        ../include/stuka/dLP/benders.h
        ../include/stuka/dLP/benders_cut.h
        ../include/stuka/dLP/benders_subproblem.h
        ../include/stuka/util/functions.h
        ../include/stuka/util/solver_factory.h
        ../include/stuka/util/timer.h
        ../include/stuka/util/callback/base_callback.h
        ../include/stuka/util/callback/function.h
        ../include/stuka/stuka.h
        ../include/stuka/constants.h
        ../include/stuka/base_solver.h
        ../include/stuka/optimize_state.h
        ../include/stuka/options.h
        LP/slack_lp.cpp
        LP/mehrotra_pc.cpp
        dLP/cre.cpp
        dLP/cre_subproblem.cpp
        dLP/critical_region.cpp
        dLP/benders_subproblem.cpp
        dLP/benders.cpp
        util/solver_factory.cpp
        util/timer.cpp
        util/functions.cpp
        util/callback/function.cpp)

add_library(stuka SHARED ${STUKA_SOURCES} ${GUROBI_SOURCES} ${HDF5_SOURCES})
target_link_libraries(stuka PRIVATE ${GUROBI_LIBRARIES} ${SPQR_LIBRARIES} ${HDF5_LIBRARIES})

set(STUKA_EXAMPLE_SOURCES
        example/simple_lp.h
        example/simple_lp.cpp
        example/simple_qp.h
        example/simple_qp.cpp
        example/borrelli_729.h
        example/borrelli_729.cc
        example/borrelli_731.h
        example/borrelli_731.cc)

add_executable(tester test.cc ${STUKA_EXAMPLE_SOURCES})
target_link_libraries(tester PRIVATE stuka)