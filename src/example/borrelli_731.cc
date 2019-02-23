//
// Created by Avinash Madavan on 7/16/18.
//

#include "borrelli_731.h"

stuka::dLP::DecomposedLinearProgram stuka::example::Borrelli731::gen() {
    typedef Eigen::Triplet<double> T;

    std::shared_ptr<Eigen::VectorXd> c0 = std::make_shared<Eigen::VectorXd>(2);
    std::shared_ptr<Eigen::VectorXd> lb0 = std::make_shared<Eigen::VectorXd>(2);
    std::shared_ptr<Eigen::VectorXd> ub0 = std::make_shared<Eigen::VectorXd>(2);

    *c0 << 0, 0;
    *lb0 << -10, -10;
    *ub0 << 10, 10;

    std::shared_ptr<Eigen::VectorXd> c1 = std::make_shared<Eigen::VectorXd>(2);
    std::shared_ptr<Eigen::SparseMatrix<double>> A_ub1 = std::make_shared<Eigen::SparseMatrix<double>>(3, 2);
    std::shared_ptr<Eigen::SparseMatrix<double>> C_ub1 = std::make_shared<Eigen::SparseMatrix<double>>(3, 2);
    std::shared_ptr<Eigen::VectorXd> b_ub1 = std::make_shared<Eigen::VectorXd>(3);
    std::shared_ptr<Eigen::VectorXd> lb1 = std::make_shared<Eigen::VectorXd>(2);
    std::shared_ptr<Eigen::VectorXd> ub1 = std::make_shared<Eigen::VectorXd>(2);

    *c1 << -2, -1;

    std::vector<T> A_ub1_triplet(6);
    A_ub1_triplet.emplace_back(T(0, 0, 1.));
    A_ub1_triplet.emplace_back(T(1, 0, 2.));
    A_ub1_triplet.emplace_back(T(2, 0, 1.));
    A_ub1_triplet.emplace_back(T(0, 1, 3.));
    A_ub1_triplet.emplace_back(T(1, 1, 1.));
    A_ub1_triplet.emplace_back(T(2, 1, 0.));
    A_ub1->setFromTriplets(A_ub1_triplet.begin(), A_ub1_triplet.end());

    std::vector<T> C_ub1_triplet(12);
    C_ub1_triplet.emplace_back(T(0, 0, 2.));
    C_ub1_triplet.emplace_back(T(1, 0, -1.));
    C_ub1_triplet.emplace_back(T(2, 0, -1.));
    C_ub1_triplet.emplace_back(T(0, 1, -1.));
    C_ub1_triplet.emplace_back(T(1, 1, 2.));
    C_ub1_triplet.emplace_back(T(2, 1, -1.));
    C_ub1->setFromTriplets(C_ub1_triplet.begin(), C_ub1_triplet.end());

    *b_ub1 << 9, 8, 4;
    *lb1 << 0, 0;
    *ub1 << INF, INF;

    std::vector<std::shared_ptr<Eigen::VectorXd>> c(2);
    std::vector<std::shared_ptr<Eigen::SparseMatrix<double>>> A_ub(2);
    std::vector<std::shared_ptr<Eigen::SparseMatrix<double>>> C_ub(1);
    std::vector<std::shared_ptr<Eigen::VectorXd>> b_ub(2);
    std::vector<std::shared_ptr<Eigen::SparseMatrix<double>>> A_eq(2);
    std::vector<std::shared_ptr<Eigen::SparseMatrix<double>>> C_eq(1);
    std::vector<std::shared_ptr<Eigen::VectorXd>> b_eq(2);
    std::vector<std::shared_ptr<Eigen::VectorXd>> lb(2);
    std::vector<std::shared_ptr<Eigen::VectorXd>> ub(2);

    c[1] = c0;
    A_ub[1] = nullptr;
    b_ub[1] = nullptr;
    A_eq[1] = nullptr;
    b_eq[1] = nullptr;
    lb[1] = lb0;
    ub[1] = ub0;

    c[0] = c1;
    A_ub[0] = A_ub1;
    C_ub[0] = C_ub1;
    b_ub[0] = b_ub1;
    A_eq[0] = nullptr;
    C_eq[0] = nullptr;
    b_eq[0] = nullptr;
    lb[0] = lb1;
    ub[0] = ub1;

    return dLP::DecomposedLinearProgram(c, A_ub, b_ub, C_ub, A_eq, b_eq, C_eq, lb, ub);
}

stuka::LP::LinearProgram stuka::example::Borrelli731::full() {
    typedef Eigen::Triplet<double> T;

    std::shared_ptr<Eigen::VectorXd> c = std::make_shared<Eigen::VectorXd>(4);
    std::shared_ptr<Eigen::VectorXd> lb = std::make_shared<Eigen::VectorXd>(4);
    std::shared_ptr<Eigen::VectorXd> ub = std::make_shared<Eigen::VectorXd>(4);

    *c << -2, -1, 0, 0;
    *lb << 0, 0, -10, -10;
    *ub << INF, INF, 10, 10;

    std::shared_ptr<Eigen::VectorXd> b_ub = std::make_shared<Eigen::VectorXd>(3);
    *b_ub << 9, 8, 4;

    std::shared_ptr<Eigen::SparseMatrix<double>> A_ub1 = std::make_shared<Eigen::SparseMatrix<double>>(3, 2);
    std::shared_ptr<Eigen::SparseMatrix<double>> C_ub1 = std::make_shared<Eigen::SparseMatrix<double>>(3, 2);

    std::vector<T> A_ub1_triplet(6);
    A_ub1_triplet.emplace_back(T(0, 0, 1.));
    A_ub1_triplet.emplace_back(T(1, 0, 2.));
    A_ub1_triplet.emplace_back(T(2, 0, 1.));
    A_ub1_triplet.emplace_back(T(0, 1, 3.));
    A_ub1_triplet.emplace_back(T(1, 1, 1.));
    A_ub1_triplet.emplace_back(T(2, 1, 0.));
    A_ub1->setFromTriplets(A_ub1_triplet.begin(), A_ub1_triplet.end());

    std::vector<T> C_ub1_triplet(12);
    C_ub1_triplet.emplace_back(T(0, 0, 2.));
    C_ub1_triplet.emplace_back(T(1, 0, -1.));
    C_ub1_triplet.emplace_back(T(2, 0, -1.));
    C_ub1_triplet.emplace_back(T(0, 1, -1.));
    C_ub1_triplet.emplace_back(T(1, 1, 2.));
    C_ub1_triplet.emplace_back(T(2, 1, -1.));
    C_ub1->setFromTriplets(C_ub1_triplet.begin(), C_ub1_triplet.end());

    std::shared_ptr<Eigen::SparseMatrix<double>> A_ub = std::make_shared<Eigen::SparseMatrix<double>>(3,4);

    A_ub->setZero();
    for (int i = 0; i < 2; ++i) {
        A_ub->startVec(i);
        for(Eigen::SparseMatrix<double>::InnerIterator it(*A_ub1, i); it; ++it)
            A_ub->insertBack(it.row(), i) = it.value();
    }
    for (int i = 0; i < 2; ++i) {
        A_ub->startVec(i+2);
        for(Eigen::SparseMatrix<double>::InnerIterator it(*C_ub1, i); it; ++it)
            A_ub->insertBack(it.row(), i+2) = it.value();
    }
    A_ub->finalize();

    LP::LinearProgram lp;
    lp.c = c;
    lp.A_ub = A_ub;
    lp.b_ub = b_ub;
    lp.lb = lb;
    lp.ub = ub;

    return lp;
}
