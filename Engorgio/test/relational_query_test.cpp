#include "ENGORGIO/comp.h"
#include "ENGORGIO/sort.h"
#include "ENGORGIO/utils.h"
#include <chrono>
#include <omp.h>
#include "math/chebyshev.h"

using namespace openfhe;
using namespace lbcrypto;

/***
 * TPC-H Query 1
    select
        l_returnflag,
        l_linestatus,
        sum(l_quantity) as sum_qty,
        sum(l_extendedprice) as sum_base_price,
        sum(l_extendedprice * (1 - l_discount)) as sum_disc_price,
        sum(l_extendedprice * (1 - l_discount) * (1 + l_tax)) as sum_charge,
        avg(l_quantity) as avg_qty,
        avg(l_extendedprice) as avg_price,
        avg(l_discount) as avg_disc,
        count(*) as count_order
    from
        lineitem
    where
        l_shipdate <= date '1998-12-01' - interval ':1' day (3)
    group by
        l_returnflag,
        l_linestatus
    order by
        l_returnflag,
        l_linestatus;
*/

double Eval_E2E_q1(int records_num)
{

    std::cout << "Relational SQL Query1 Test: " << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Records: " << records_num << std::endl;
    CCParams<CryptoContextCKKSRNS> parameters;

    parameters.SetSecurityLevel(HEStd_128_classic);
#if NATIVEINT == 128
    usint scalingModSize = 78;
    usint firstModSize = 89;
#else
    usint scalingModSize = 50;
    usint firstModSize = 60;
#endif
    SecretKeyDist secretKeyDist = SPARSE_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);
    parameters.SetScalingModSize(scalingModSize);
    parameters.SetFirstModSize(firstModSize);
    std::uint32_t polyDegree = 119;

    std::uint32_t multDepth = 27;
    parameters.SetMultiplicativeDepth(multDepth);
    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);
    std::uniform_int_distribution<std::uint32_t> shipdate_message(0, 10592 + 100);
    std::uniform_int_distribution<uint64_t> quantity_message(0, 8);
    std::uniform_real_distribution<double> extendedprice_message(0., 8.);
    std::uniform_real_distribution<double> discount_message(0., 1.);
    std::uniform_real_distribution<double> tax_message(0., 1.);
    std::uniform_int_distribution<std::uint32_t> bianry_message(0, 1);
    double precision = (1 << (8 - 1)) - 1;
    std::cout << precision << std::endl;
    double lowerBound = -precision;
    double upperBound = precision;
    double bound = 3;
    auto keyPair = cc->KeyGen();
    usint ringDim = cc->GetRingDimension();
    int length = ringDim / 2;
    int Bg = 7;
    cc->EvalMultKeyGen(keyPair.secretKey);

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<int> message(0, precision);
    std::vector<int> rotstep;
    for (int i = 1; i < length; i *= 2)
    {
        rotstep.push_back(i);
        rotstep.push_back(-i);
    }

    cc->EvalRotateKeyGen(keyPair.secretKey, rotstep);

    std::vector<double> ship_date_data(length), ship_date_predict(length, 10592), returnflag(length), linestatus(length), returnflag_predicate_Y(length, 100),
        returnflag_predicate_N(length, 0), linestatus_predicate_Y(length, 100), linestatus_predicate_N(length, 0),
        quantity(length), extendedprice(length), discount(length), tax(length);

    std::vector<Ciphertext<lbcrypto::DCRTPoly>> ship_data_ciphers, ship_predict_ciphers;
    Ciphertext<lbcrypto::DCRTPoly> returnflag_ciphers, linestatus_ciphers, returnflag_predicate_Y_ciphers,
        returnflag_predicate_N_ciphers, linestatus_predicate_Y_ciphers, linestatus_predicate_N_ciphers, qty_ciphers;

    for (int i = 0; i < length; i++)
    {
        ship_date_data[i] = double(shipdate_message(engine));
        returnflag[i] = bianry_message(engine) * 100;
        linestatus[i] = bianry_message(engine) * 100;
        quantity[i] = quantity_message(engine);
        extendedprice[i] = extendedprice_message(engine);
        discount[i] = discount_message(engine);
        tax[i] = tax_message(engine);
    }
    auto ship_date_data_quant = quantization(ship_date_data, 16, Bg);
    auto ship_date_predict_quant = quantization(ship_date_predict, 16, Bg);

    int block = ship_date_data_quant.size();
    size_t encodedLength = ship_date_data.size();
    for (int i = 0; i < block; i++)
    {
        Plaintext plain_a = cc->MakeCKKSPackedPlaintext(ship_date_data_quant[i]);
        Plaintext plain_b = cc->MakeCKKSPackedPlaintext(ship_date_predict_quant[i]);
        ship_data_ciphers.push_back(cc->Encrypt(keyPair.publicKey, plain_a));
        ship_predict_ciphers.push_back(cc->Encrypt(keyPair.publicKey, plain_b));
    }
    Plaintext plain_linestatus, plain_returnflag, plain_returnflag_Y, plain_returnflag_N, plain_linestatus_Y, plain_linestatus_N, plain_qty;

    plain_qty = cc->MakeCKKSPackedPlaintext(quantity);
    qty_ciphers = cc->Encrypt(keyPair.publicKey, plain_qty);

    plain_linestatus = cc->MakeCKKSPackedPlaintext(linestatus);
    plain_linestatus_Y = cc->MakeCKKSPackedPlaintext(linestatus_predicate_Y);
    plain_linestatus_N = cc->MakeCKKSPackedPlaintext(linestatus_predicate_N);

    plain_returnflag = cc->MakeCKKSPackedPlaintext(returnflag);
    plain_returnflag_Y = cc->MakeCKKSPackedPlaintext(returnflag_predicate_Y);
    plain_returnflag_N = cc->MakeCKKSPackedPlaintext(returnflag_predicate_N);

    returnflag_ciphers = cc->Encrypt(keyPair.publicKey, plain_returnflag);
    linestatus_ciphers = cc->Encrypt(keyPair.publicKey, plain_linestatus);
    returnflag_predicate_Y_ciphers = cc->Encrypt(keyPair.publicKey, plain_linestatus_Y);
    linestatus_predicate_Y_ciphers = cc->Encrypt(keyPair.publicKey, plain_returnflag_Y);
    linestatus_predicate_N_ciphers = cc->Encrypt(keyPair.publicKey, plain_returnflag_N);

    Ciphertext<lbcrypto::DCRTPoly> filter_res_YY, filter_res_YN, filter_res_NY, filter_res_NN;
    double filtering_time = 0,
           aggregation_time = 0;
    std::chrono::system_clock::time_point start, end;
    Ciphertext<lbcrypto::DCRTPoly> comp_res_ge, comp_res_eq, comp_returnflag_YY, comp_returnflag_YN, comp_linestatus_YY, comp_linestatus_YN;
    // filtering (can be parallel)
    start = std::chrono::system_clock::now();

    comp_greater_than_modular(ship_predict_ciphers, ship_data_ciphers, precision, polyDegree, comp_res_ge, keyPair.secretKey);

    comp_equal(returnflag_ciphers, returnflag_predicate_Y_ciphers, precision, polyDegree, comp_returnflag_YY);

    comp_equal(linestatus_ciphers, linestatus_predicate_Y_ciphers, precision, polyDegree, comp_linestatus_YY);

    auto comp_returnflag_YY_neg = cc->EvalNegate(comp_returnflag_YY);
    comp_returnflag_YN = cc->EvalAdd(comp_returnflag_YY_neg, 1.0);
    auto comp_linestatus_YY_neg = cc->EvalNegate(comp_linestatus_YY);
    comp_linestatus_YN = cc->EvalAdd(comp_linestatus_YY_neg, 1.0);

    filter_res_YY = cc->EvalMult(comp_res_ge, comp_returnflag_YY);
    filter_res_YY = cc->EvalMult(filter_res_YY, comp_linestatus_YY);

    filter_res_YN = cc->EvalMult(comp_res_ge, comp_returnflag_YY);
    filter_res_YN = cc->EvalMult(filter_res_YN, comp_linestatus_YN);

    filter_res_NY = cc->EvalMult(comp_res_ge, comp_returnflag_YN);
    filter_res_NY = cc->EvalMult(filter_res_NY, comp_linestatus_YY);

    filter_res_NN = cc->EvalMult(comp_res_ge, comp_returnflag_YN);
    filter_res_NN = cc->EvalMult(filter_res_NN, comp_linestatus_YN);

    end = std::chrono::system_clock::now();
    filtering_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() * records_num / length * 3;
    printf("filter time = %f ms\n", filtering_time);

    std::vector<uint64_t> plain_filter_res_YY(length, 0), plain_filter_res_YN(length, 0), plain_filter_res_NY(length, 0), plain_filter_res_NN(length, 0);
    uint64_t plain_agg_res_YY = 0, plain_agg_res_YN = 0, plain_agg_res_NY = 0, plain_agg_res_NN = 0;

    for (size_t i = 0; i < length; i++)
    {
        if (ship_date_data[i] < 10592)
        {
            if (returnflag[i] == 100)
            {
                if (linestatus[i] == 100)
                {
                    plain_filter_res_YY[i] = 1;
                    plain_agg_res_YY += quantity[i];
                }
                else
                {
                    plain_filter_res_YN[i] = 1;
                    plain_agg_res_YN += quantity[i];
                }
            }
            else
            {
                if (linestatus[i] == 100)
                {
                    plain_filter_res_NY[i] = 1;
                    plain_agg_res_NY += quantity[i];
                }
                else
                {
                    plain_filter_res_NN[i] = 1;
                    plain_agg_res_NN += quantity[i];
                }
            }
        }
    }

    std::cout << "Filtering finish" << std::endl;

    std::cout << "Aggregation :" << std::endl;
    std::cout << "Aggregating quanlity, Taking SUM(quantlity) as an example.." << std::endl;

    Ciphertext<lbcrypto::DCRTPoly> sum_qty_cipher_YY, sum_qty_cipher_YN, sum_qty_cipher_NY, sum_qty_cipher_NN;
    sum_qty_cipher_YY = cc->EvalMult(qty_ciphers, filter_res_YY);
    sum_qty_cipher_YN = cc->EvalMult(qty_ciphers, filter_res_YN);
    sum_qty_cipher_NY = cc->EvalMult(qty_ciphers, filter_res_NY);
    sum_qty_cipher_NN = cc->EvalMult(qty_ciphers, filter_res_NN);

    start = std::chrono::system_clock::now();
    int logrow = log2(length);
    Ciphertext<lbcrypto::DCRTPoly> temp, rot_temp;
    for (size_t i = 0; i < logrow; i++)
    {
        int step = 1 << (logrow - i - 1);
        temp = sum_qty_cipher_YY;
        rot_temp = cc->EvalRotate(temp, step);
        sum_qty_cipher_YY = cc->EvalAdd(sum_qty_cipher_YY, rot_temp);

        temp = sum_qty_cipher_YN;
        rot_temp = cc->EvalRotate(temp, step);
        sum_qty_cipher_YN = cc->EvalAdd(sum_qty_cipher_YN, rot_temp);

        temp = sum_qty_cipher_NY;
        rot_temp = cc->EvalRotate(temp, step);
        sum_qty_cipher_NY = cc->EvalAdd(sum_qty_cipher_NY, rot_temp);

        temp = sum_qty_cipher_NN;
        rot_temp = cc->EvalRotate(temp, step);
        sum_qty_cipher_NN = cc->EvalAdd(sum_qty_cipher_NN, rot_temp);
    }
    end = std::chrono::system_clock::now();

    double ta = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000;
    double agg_time = ta * records_num / length;
    printf("Agg time = %f ms\n", ta * records_num / length);

    Plaintext plaintextDec_agg_resultYY, plaintextDec_agg_resultYN, plaintextDec_agg_resultNY, plaintextDec_agg_resultNN;
    cc->Decrypt(keyPair.secretKey, sum_qty_cipher_YY, &plaintextDec_agg_resultYY);
    cc->Decrypt(keyPair.secretKey, sum_qty_cipher_YN, &plaintextDec_agg_resultYN);
    cc->Decrypt(keyPair.secretKey, sum_qty_cipher_NY, &plaintextDec_agg_resultNY);
    cc->Decrypt(keyPair.secretKey, sum_qty_cipher_NN, &plaintextDec_agg_resultNN);

    plaintextDec_agg_resultYY->SetLength(encodedLength);
    plaintextDec_agg_resultYN->SetLength(encodedLength);
    plaintextDec_agg_resultNY->SetLength(encodedLength);
    plaintextDec_agg_resultNN->SetLength(encodedLength);

    std::vector<std::complex<double>> agg_resultYY = plaintextDec_agg_resultYY->GetCKKSPackedValue();
    std::vector<std::complex<double>> agg_resultYN = plaintextDec_agg_resultYN->GetCKKSPackedValue();
    std::vector<std::complex<double>> agg_resultNY = plaintextDec_agg_resultNY->GetCKKSPackedValue();
    std::vector<std::complex<double>> agg_resultNN = plaintextDec_agg_resultNN->GetCKKSPackedValue();

    std::cout << "Query Evaluation Time: " << filtering_time + agg_time << " ms" << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;

    std::cout << "Encrypted query result: " << std::endl;
    std::cout << std::setw(16) << "returnfalg" << "|" << std::setw(16) << "linestatus" << "|" << std::setw(16) << "sum_qty" << std::endl;
    std::cout << std::setw(16) << "Y" << "|" << std::setw(16) << "Y" << "|" << std::setw(16) << std::round(agg_resultYY[0].real()) << std::endl;
    std::cout << std::setw(16) << "Y" << "|" << std::setw(16) << "N" << "|" << std::setw(16) << std::round(agg_resultYN[0].real()) << std::endl;
    std::cout << std::setw(16) << "N" << "|" << std::setw(16) << "Y" << "|" << std::setw(16) << std::round(agg_resultNY[0].real()) << std::endl;
    std::cout << std::setw(16) << "N" << "|" << std::setw(16) << "N" << "|" << std::setw(16) << std::round(agg_resultNN[0].real()) << std::endl;

    std::cout << "Plain query result: " << std::endl;
    std::cout << std::setw(16) << "returnflag" << "|" << std::setw(16) << "linestatus" << "|" << std::setw(16) << "sum_qty" << std::endl;
    std::cout << std::setw(16) << "Y" << "|" << std::setw(16) << "Y" << "|" << std::setw(16) << plain_agg_res_YY << std::endl;
    std::cout << std::setw(16) << "Y" << "|" << std::setw(16) << "N" << "|" << std::setw(16) << plain_agg_res_YN << std::endl;
    std::cout << std::setw(16) << "N" << "|" << std::setw(16) << "Y" << "|" << std::setw(16) << plain_agg_res_NY << std::endl;
    std::cout << std::setw(16) << "N" << "|" << std::setw(16) << "N" << "|" << std::setw(16) << plain_agg_res_NN << std::endl;

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    return filtering_time + agg_time;
}

/*
    TPC-H Query 12
    select
        l_shipmode,
        sum(case
            when o_orderpriority = '1-URGENT'
                or o_orderpriority = '2-HIGH'
                then 1
            else 0
        end) as high_line_count,
        sum(case
            when o_orderpriority <> '1-URGENT'
                and o_orderpriority <> '2-HIGH'
                then 1
            else 0
        end) as low_line_count
    from
        orders,
        lineitem
    where
        o_orderkey = l_orderkey
        and l_shipmode in (':1', ':2')
        and l_commitdate < l_receiptdate
        and l_shipdate < l_commitdate
        and l_receiptdate >= date ':3'
        and l_receiptdate < date ':3' + interval '1' year
    group by
        l_shipmode
    order by
        l_shipmode;
    Consider the joined table
*/

double Eval_E2E_q12(int records_num)
{

    std::cout << "Relational SQL Query12 Test: " << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Records: " << records_num << std::endl;
    CCParams<CryptoContextCKKSRNS> parameters;

    parameters.SetSecurityLevel(HEStd_128_classic);
#if NATIVEINT == 128
    usint scalingModSize = 78;
    usint firstModSize = 89;
#else
    usint scalingModSize = 50;
    usint firstModSize = 60;
#endif
    SecretKeyDist secretKeyDist = SPARSE_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);
    parameters.SetScalingModSize(scalingModSize);
    parameters.SetFirstModSize(firstModSize);
    std::uint32_t polyDegree = 119;

    std::uint32_t multDepth = 30;
    parameters.SetMultiplicativeDepth(multDepth);
    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);
    std::uniform_int_distribution<std::uint32_t> shipmode_message(1, 10);
    std::uniform_int_distribution<std::uint32_t> shipdate_message(0, 15000);
    std::uniform_int_distribution<std::uint32_t> receiptdate_message(0, 15000);
    std::uniform_int_distribution<std::uint32_t> commitdate_message(0, 15000);
    // orderpriority \in ('1-URGENT', '2-HIGH', '3-MEDIUM', '4-NOT SPECIFIED', '5-LOW')
    std::uniform_int_distribution<uint64_t> orderpriority_message(1, 5);

    double precision = (1 << (8 - 1)) - 1;
    std::cout << precision << std::endl;
    double lowerBound = -precision;
    double upperBound = precision;
    double bound = 3;
    auto keyPair = cc->KeyGen();
    usint ringDim = cc->GetRingDimension();
    int length = ringDim / 2;
    int Bg = 7;
    cc->EvalMultKeyGen(keyPair.secretKey);

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<int> message(0, precision);
    std::vector<int> rotstep;
    for (int i = 1; i < length; i *= 2)
    {
        rotstep.push_back(i);
        rotstep.push_back(-i);
    }

    cc->EvalRotateKeyGen(keyPair.secretKey, rotstep);

    std::vector<double> shipdate(length), commitdate(length), receiptdate(length), orderpriority(length), shipmode(length),
        returnflag_predicate_N(length, 0), linestatus_predicate_Y(length, 100), linestatus_predicate_N(length, 0),
        quantity(length), extendedprice(length), discount(length), tax(length), predicate_date1(length, 10000),
        predicate_date2(length, 13000), predicate_mail(length, 10), predicate_ship(length, 20), predicate_urgent(length, 10), predicate_high(length, 20);

    std::vector<Ciphertext<lbcrypto::DCRTPoly>> shipdate_ciphers, commitdate_ciphers, receiptdate_ciphers,
        predicate_date_cipher1, predicate_date_cipher2;
    Ciphertext<lbcrypto::DCRTPoly> orderpriority_ciphers, shipmode_ciphers, predicate_mail_cipher, predicate_ship_cipher,
        predicate_urgent_cipher, predicate_high_cipher, filter_res_mail, filter_res_ship, order_res, res_mail_order,
        res_ship_order, count_mail, count_ship, count_mail_order, count_ship_order, agg_mail, agg_ship, agg_mail_order, agg_ship_order;

    for (int i = 0; i < length; i++)
    {
        shipdate[i] = shipdate_message(engine);
        commitdate[i] = commitdate_message(engine);
        receiptdate[i] = receiptdate_message(engine);
        shipmode[i] = shipmode_message(engine) * 10;
        orderpriority[i] = orderpriority_message(engine) * 10;
    }
    auto shipdate_data_quant = quantization(shipdate, 16, Bg);
    auto commitdate_quant = quantization(commitdate, 16, Bg);
    auto receiptdate_quant = quantization(receiptdate, 16, Bg);
    auto predicate_date1_quant = quantization(predicate_date1, 16, Bg);
    auto predicate_date2_quant = quantization(predicate_date2, 16, Bg);

    int block = shipdate_data_quant.size();
    size_t encodedLength = shipdate.size();
    for (int i = 0; i < block; i++)
    {
        Plaintext plain_a = cc->MakeCKKSPackedPlaintext(shipdate_data_quant[i]);
        Plaintext plain_b = cc->MakeCKKSPackedPlaintext(commitdate_quant[i]);
        Plaintext plain_c = cc->MakeCKKSPackedPlaintext(receiptdate_quant[i]);
        Plaintext plain_d = cc->MakeCKKSPackedPlaintext(predicate_date1_quant[i]);
        Plaintext plain_e = cc->MakeCKKSPackedPlaintext(predicate_date2_quant[i]);
        shipdate_ciphers.push_back(cc->Encrypt(keyPair.publicKey, plain_a));
        commitdate_ciphers.push_back(cc->Encrypt(keyPair.publicKey, plain_b));
        receiptdate_ciphers.push_back(cc->Encrypt(keyPair.publicKey, plain_c));
        predicate_date_cipher1.push_back(cc->Encrypt(keyPair.publicKey, plain_d));
        predicate_date_cipher2.push_back(cc->Encrypt(keyPair.publicKey, plain_e));
    }
    Plaintext plain_shipmode, plain_orderpriority, plain_predicate_mail, plain_predicate_ship, plain_predicate_urgent, plain_predicate_high;

    plain_shipmode = cc->MakeCKKSPackedPlaintext(shipmode);
    shipmode_ciphers = cc->Encrypt(keyPair.publicKey, plain_shipmode);

    plain_orderpriority = cc->MakeCKKSPackedPlaintext(orderpriority);
    orderpriority_ciphers = cc->Encrypt(keyPair.publicKey, plain_orderpriority);

    plain_predicate_mail = cc->MakeCKKSPackedPlaintext(predicate_mail);
    predicate_mail_cipher = cc->Encrypt(keyPair.publicKey, plain_predicate_mail);

    plain_predicate_ship = cc->MakeCKKSPackedPlaintext(predicate_ship);
    predicate_ship_cipher = cc->Encrypt(keyPair.publicKey, plain_predicate_ship);

    plain_predicate_urgent = cc->MakeCKKSPackedPlaintext(predicate_urgent);
    predicate_urgent_cipher = cc->Encrypt(keyPair.publicKey, plain_predicate_urgent);

    plain_predicate_high = cc->MakeCKKSPackedPlaintext(predicate_high);
    predicate_high_cipher = cc->Encrypt(keyPair.publicKey, plain_predicate_high);

    double filtering_time = 0;
    std::chrono::system_clock::time_point start, end;
    Ciphertext<lbcrypto::DCRTPoly> pre_res;
    // filtering (can be parallel)
    start = std::chrono::system_clock::now();
    comp_greater_than_modular(receiptdate_ciphers, commitdate_ciphers, precision, polyDegree, filter_res_mail, keyPair.secretKey);

    comp_greater_than_modular(commitdate_ciphers, shipdate_ciphers, precision, polyDegree, pre_res, keyPair.secretKey);

    filter_res_mail = cc->EvalMult(filter_res_mail, pre_res);

    comp_greater_than_modular(receiptdate_ciphers, predicate_date_cipher1, precision, polyDegree, pre_res, keyPair.secretKey);

    filter_res_mail = cc->EvalMult(filter_res_mail, pre_res);

    comp_greater_than_modular(predicate_date_cipher2, receiptdate_ciphers, precision, polyDegree, pre_res, keyPair.secretKey);

    filter_res_mail = cc->EvalMult(filter_res_mail, pre_res);

    filter_res_ship = filter_res_mail;

    comp_equal(shipmode_ciphers, predicate_mail_cipher, precision, polyDegree, pre_res);

    filter_res_mail = cc->EvalMult(filter_res_mail, pre_res);

    comp_equal(shipmode_ciphers, predicate_ship_cipher, precision, polyDegree, pre_res);

    filter_res_ship = cc->EvalMult(filter_res_ship, pre_res);

    comp_equal(orderpriority_ciphers, predicate_urgent_cipher, precision, polyDegree, order_res);

    comp_equal(orderpriority_ciphers, predicate_high_cipher, precision, polyDegree, pre_res);

    auto neg_pre_res = cc->EvalNegate(pre_res);
    auto neg_order_res = cc->EvalNegate(order_res);
    neg_pre_res = cc->EvalAdd(neg_pre_res, 1.0);
    neg_order_res = cc->EvalAdd(neg_order_res, 1.0);

    auto YY_pre_order = cc->EvalMult(pre_res, order_res);
    auto NY_pre_order = cc->EvalMult(pre_res, neg_order_res);
    auto YN_pre_order = cc->EvalMult(neg_pre_res, order_res);
    order_res = cc->EvalAdd(YY_pre_order, NY_pre_order);
    order_res = cc->EvalAdd(order_res, YN_pre_order);

    count_mail_order = cc->EvalMult(filter_res_mail, order_res);

    count_ship_order = cc->EvalMult(filter_res_ship, order_res);

    count_mail = cc->EvalMult(filter_res_mail, filter_res_mail);

    count_ship = cc->EvalMult(filter_res_ship, filter_res_ship);

    int logrow = log2(length);

    Ciphertext<lbcrypto::DCRTPoly> temp, rot_temp;

    agg_mail = count_mail;
    agg_ship = count_ship;
    agg_mail_order = count_mail_order;
    agg_ship_order = count_ship_order;
    for (size_t i = 0; i < logrow; i++)
    {
        int step = 1 << (logrow - i - 1);
        temp = agg_mail;
        rot_temp = cc->EvalRotate(temp, step);
        agg_mail = cc->EvalAdd(agg_mail, rot_temp);

        temp = agg_ship;
        rot_temp = cc->EvalRotate(temp, step);
        agg_ship = cc->EvalAdd(agg_ship, rot_temp);

        temp = agg_mail_order;
        rot_temp = cc->EvalRotate(temp, step);
        agg_mail_order = cc->EvalAdd(agg_mail_order, rot_temp);

        temp = agg_ship_order;
        rot_temp = cc->EvalRotate(temp, step);
        agg_ship_order = cc->EvalAdd(agg_ship_order, rot_temp);
    }

    end = std::chrono::system_clock::now();

    filtering_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() * records_num / length;
    printf("filter time = %f ms\n", filtering_time);
    std::cout << "Filtering finish" << std::endl;
    std::vector<uint64_t> plain_filter_res_mail(length, 0), plain_filter_res_ship(length, 0), plain_filter_order(length, 0);
    std::vector<uint64_t> plain_res_mail_order(length, 0), plain_res_ship_order(length, 0);
    uint64_t agg_mail_res = 0, agg_mail_order_res = 0, agg_ship_res = 0, agg_ship_order_res = 0;
    bool ress;
    start = std::chrono::system_clock::now();
    for (size_t i = 0; i < length; i++)
    {
        if (commitdate[i] < receiptdate[i] && shipdate[i] < commitdate[i] && receiptdate[i] > predicate_date1[i] && receiptdate[i] < predicate_date2[i])
        {
            ress = true;
        }
        else
        {
            ress = false;
        }

        if (orderpriority[i] == 10 || orderpriority[i] == 20)
        {
            plain_filter_order[i] = 1;
        }

        if (ress && shipmode[i] == 10)
        {
            plain_filter_res_mail[i] = 1;
            agg_mail_res += 1;
            if (plain_filter_order[i] == 1)
            {
                plain_res_mail_order[i] = 1;
                agg_mail_order_res += 1;
            }
        }
        if (ress && shipmode[i] == 20)
        {
            plain_filter_res_ship[i] = 1;
            agg_ship_res += 1;
            if (plain_filter_order[i] == 1)
            {
                plain_res_ship_order[i] = 1;
                agg_ship_order_res += 1;
            }
        }
    }
    end = std::chrono::system_clock::now();

    double ta = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000;
    double agg_time = ta * records_num / length;
    printf("Agg time = %f ms\n", agg_time);

    Plaintext query_res_mail, query_res_ship, query_res_mail_order, query_res_ship_order;
    cc->Decrypt(keyPair.secretKey, agg_mail, &query_res_mail);
    cc->Decrypt(keyPair.secretKey, agg_ship, &query_res_ship);
    cc->Decrypt(keyPair.secretKey, agg_mail_order, &query_res_mail_order);
    cc->Decrypt(keyPair.secretKey, agg_ship_order, &query_res_ship_order);

    query_res_mail->SetLength(encodedLength);
    query_res_ship->SetLength(encodedLength);
    query_res_mail_order->SetLength(encodedLength);
    query_res_ship_order->SetLength(encodedLength);

    std::vector<std::complex<double>> res_mail_dec = query_res_mail->GetCKKSPackedValue();
    std::vector<std::complex<double>> res_ship_dec = query_res_ship->GetCKKSPackedValue();
    std::vector<std::complex<double>> res_mail_order_dec = query_res_mail_order->GetCKKSPackedValue();
    std::vector<std::complex<double>> res_ship_order_dec = query_res_ship_order->GetCKKSPackedValue();

    std::cout << "Query Evaluation Time: " << filtering_time + agg_time << " ms" << std::endl;
    std::cout << "Encrypted result: " << std::endl;
    std::cout << std::setw(12) << "shipmode" << "|" << std::setw(16) << "high_line_count" << "|" << std::setw(16) << "low_line_count" << std::endl;
    std::cout << std::setw(12) << "MAIL" << "|" << std::setw(16) << res_mail_order_dec[0].real() << "|" << std::setw(16) << res_mail_dec[0].real() - res_mail_order_dec[0].real() << std::endl;
    std::cout << std::setw(12) << "SHIP" << "|" << std::setw(16) << res_ship_order_dec[0].real() << "|" << std::setw(16) << res_ship_dec[0].real() - res_ship_order_dec[0].real() << std::endl;

    std::cout << "Plain result: " << std::endl;
    std::cout << std::setw(12) << "shipmode" << "|" << std::setw(16) << "high_line_count" << "|" << std::setw(16) << "low_line_count" << std::endl;
    std::cout << std::setw(12) << "MAIL" << "|" << std::setw(16) << agg_mail_order_res << "|" << std::setw(16) << agg_mail_res - agg_mail_order_res << std::endl;
    std::cout << std::setw(12) << "SHIP" << "|" << std::setw(16) << agg_ship_order_res << "|" << std::setw(16) << agg_ship_res - agg_ship_order_res << std::endl;

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    return filtering_time + agg_time;
}

/***
 * TPC-H Query 6
 * select
        sum(l_extendedprice * l_discount) as revenue
    from
        lineitem
    where
        l_shipdate >= date ':1'
        and l_shipdate < date ':1' + interval '1' year
        and l_discount between :2 - 0.01 and :2 + 0.01
        and l_quantity < :3;

*/

// 16-bit precision
double Eval_E2E_q6(int records_num)
{

    std::cout << "Relational SQL Query6 Test: " << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Records: " << records_num << std::endl;
    CCParams<CryptoContextCKKSRNS> parameters;

    parameters.SetSecurityLevel(HEStd_128_classic);
#if NATIVEINT == 128
    usint scalingModSize = 78;
    usint firstModSize = 89;
#else
    usint scalingModSize = 50;
    usint firstModSize = 60;
#endif
    SecretKeyDist secretKeyDist = SPARSE_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);
    parameters.SetScalingModSize(scalingModSize);
    parameters.SetFirstModSize(firstModSize);
    std::uint32_t polyDegree = 119;

    std::uint32_t multDepth = 30;
    parameters.SetMultiplicativeDepth(multDepth);
    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);
    std::uniform_int_distribution<uint32_t> shipdate_message(10000, 15000);
    std::uniform_int_distribution<uint32_t> discount_message(2000, 9000);
    std::uniform_int_distribution<uint32_t> quantity_message(2400, 2500);
    std::uniform_int_distribution<uint64_t> revenue_message(0, 100);

    double precision = (1 << (8 - 1)) - 1;
    std::cout << precision << std::endl;
    double lowerBound = -precision;
    double upperBound = precision;
    double bound = 3;
    auto keyPair = cc->KeyGen();
    usint ringDim = cc->GetRingDimension();
    std::cout << "N = " << ringDim << std::endl;

    int length = ringDim / 2;
    int Bg = 7;
    cc->EvalMultKeyGen(keyPair.secretKey);

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<int> message(0, precision);
    std::vector<int> rotstep;
    for (int i = 1; i < length; i *= 2)
    {
        rotstep.push_back(i);
        rotstep.push_back(-i);
    }

    cc->EvalRotateKeyGen(keyPair.secretKey, rotstep);


    
    auto cryptoParams = cc->GetCryptoParameters();
    auto params = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cryptoParams);

    
    auto moduli = params->GetElementParams()->GetParams();

    std::cout << "Levels = " << moduli.size() << std::endl;

    double total_bits = 0;
    for (size_t i = 0; i < moduli.size(); i++) {
        auto qi = moduli[i]->GetModulus();
        size_t bits = qi.GetMSB();

        std::cout << "qi[" << i << "] bits = " << bits << std::endl;
        total_bits += bits;
    }

    std::cout << "Total Q bits ≈ " << total_bits << std::endl;
    std::cout << "Levels = " << moduli.size() << std::endl;
    
    auto q0 = moduli[0]->GetModulus();
    std::cout << "First modulus q0 = " << q0
            << " (bits = " << q0.GetMSB() << ")" << std::endl;

    
    auto q1 = moduli[1]->GetModulus();
    std::cout << "Scaling modulus qi = " << q1
            << " (bits = " << q1.GetMSB() << ")" << std::endl;


    std::vector<double> revenue(length), ship_date(length), discount(length), quantity(length), predicate1_value(length, 10592),
        predicate2_value(length, 10957), predicate3_value(length, 3500), predicate4_value(length, 4500), predicate5_value(length, 2450);

    std::vector<Ciphertext<lbcrypto::DCRTPoly>> shipdate_ciphers, discount_ciphers, quantity_ciphers, predicate_value_cipher1,
        predicate_value_cipher2,
        predicate_value_cipher3,
        predicate_value_cipher4,
        predicate_value_cipher5;
    Ciphertext<lbcrypto::DCRTPoly> revenue_cipher, filter_res, result;

    for (int i = 0; i < length; i++)
    {
        revenue[i] = revenue_message(engine);
        ship_date[i] = shipdate_message(engine);
        discount[i] = discount_message(engine);
        quantity[i] = quantity_message(engine);
    }
    auto shipdate_data_quant = quantization(ship_date, 16, Bg);
    auto discount_quant = quantization(discount, 16, Bg);
    auto quantity_quant = quantization(quantity, 16, Bg);
    auto predicate1_value_quant = quantization(predicate1_value, 16, Bg);
    auto predicate2_value_quant = quantization(predicate2_value, 16, Bg);
    auto predicate3_value_quant = quantization(predicate3_value, 16, Bg);
    auto predicate4_value_quant = quantization(predicate4_value, 16, Bg);
    auto predicate5_value_quant = quantization(predicate5_value, 16, Bg);

    int block = shipdate_data_quant.size();
    size_t encodedLength = ship_date.size();
    for (int i = 0; i < block; i++)
    {
        Plaintext plain_a = cc->MakeCKKSPackedPlaintext(shipdate_data_quant[i]);
        Plaintext plain_b = cc->MakeCKKSPackedPlaintext(discount_quant[i]);
        Plaintext plain_c = cc->MakeCKKSPackedPlaintext(quantity_quant[i]);
        Plaintext plain_d = cc->MakeCKKSPackedPlaintext(predicate1_value_quant[i]);
        Plaintext plain_e = cc->MakeCKKSPackedPlaintext(predicate2_value_quant[i]);
        Plaintext plain_f = cc->MakeCKKSPackedPlaintext(predicate3_value_quant[i]);
        Plaintext plain_g = cc->MakeCKKSPackedPlaintext(predicate4_value_quant[i]);
        Plaintext plain_h = cc->MakeCKKSPackedPlaintext(predicate5_value_quant[i]);
        shipdate_ciphers.push_back(cc->Encrypt(keyPair.publicKey, plain_a));
        discount_ciphers.push_back(cc->Encrypt(keyPair.publicKey, plain_b));
        quantity_ciphers.push_back(cc->Encrypt(keyPair.publicKey, plain_c));
        predicate_value_cipher1.push_back(cc->Encrypt(keyPair.publicKey, plain_d));
        predicate_value_cipher2.push_back(cc->Encrypt(keyPair.publicKey, plain_e));
        predicate_value_cipher3.push_back(cc->Encrypt(keyPair.publicKey, plain_f));
        predicate_value_cipher4.push_back(cc->Encrypt(keyPair.publicKey, plain_g));
        predicate_value_cipher5.push_back(cc->Encrypt(keyPair.publicKey, plain_h));
    }
    Plaintext plain_revenue;

    plain_revenue = cc->MakeCKKSPackedPlaintext(revenue);
    revenue_cipher = cc->Encrypt(keyPair.publicKey, plain_revenue);

    double filtering_time = 0, aggregation_time = 0;
    std::chrono::system_clock::time_point start, end, start0, end0;
    Ciphertext<lbcrypto::DCRTPoly> pre_res;

    // filtering (can be parallel)
    start = std::chrono::system_clock::now();

    comp_greater_than_modular(shipdate_ciphers, predicate_value_cipher1, precision, polyDegree, filter_res, keyPair.secretKey);

    comp_greater_than_modular(predicate_value_cipher2, shipdate_ciphers, precision, polyDegree, pre_res, keyPair.secretKey);

    filter_res = cc->EvalMult(pre_res, filter_res);

    comp_greater_than_modular(discount_ciphers, predicate_value_cipher3, precision, polyDegree, pre_res, keyPair.secretKey);

    filter_res = cc->EvalMult(pre_res, filter_res);

    comp_greater_than_modular(predicate_value_cipher4, discount_ciphers, precision, polyDegree, pre_res, keyPair.secretKey);

    filter_res = cc->EvalMult(pre_res, filter_res);

    comp_greater_than_modular(predicate_value_cipher5, quantity_ciphers, precision, polyDegree, pre_res, keyPair.secretKey);

    filter_res = cc->EvalMult(pre_res, filter_res);
    result = cc->EvalMult(filter_res, revenue_cipher);

    int logrow = log2(length);
    std::vector<uint64_t> plain_filter_res(length);
    uint64_t plain_agg_res = 0;
    for (size_t i = 0; i < length; i++)
    {
        if (ship_date[i] > predicate1_value[0] && ship_date[i] < predicate2_value[0] &&
            discount[i] > predicate3_value[0] && discount[i] < predicate4_value[0] && quantity[i] < predicate5_value[0])
        {
            plain_filter_res[i] = 1.;
            plain_agg_res += revenue[i];
        }
        else
        {
            plain_filter_res[i] = 0.;
        }
    }
    end = std::chrono::system_clock::now();

    filtering_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() * records_num / length;
    printf("filter time = %f ms\n", filtering_time);
    std::cout << "Filtering finish" << std::endl;
    Ciphertext<lbcrypto::DCRTPoly> temp, rot_temp;
    // aggregation
    start = std::chrono::system_clock::now();
    for (size_t i = 0; i < logrow; i++)
    {
        int step = 1 << (logrow - i - 1);
        temp = result;
        start0 = std::chrono::system_clock::now();
        rot_temp = cc->EvalRotate(temp, step);
        end0 = std::chrono::system_clock::now();
        std::cout << "Rotate time = " << std::chrono::duration_cast<std::chrono::milliseconds>(end0 - start0).count()
                  << " ms" << std::endl;
        result = cc->EvalAdd(result, rot_temp);
    }

    end = std::chrono::system_clock::now();

    aggregation_time += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() * records_num / length;
    printf("aggregation_time time = %f ms\n", aggregation_time);
    Plaintext agg_result, query_res_ship, query_res_mail_order, query_res_ship_order;
    cc->Decrypt(keyPair.secretKey, result, &agg_result);

    agg_result->SetLength(encodedLength);

    std::vector<std::complex<double>> res_dec = agg_result->GetCKKSPackedValue();

    std::cout << "Query Evaluation Time: " << filtering_time + aggregation_time << " ms" << std::endl;

    std::cout << "Encrypted query result: " << std::endl;
    std::cout << std::setw(12) << "revenue" << std::endl;
    std::cout << std::setw(12) << std::round(res_dec[0].real()) << std::endl;
    std::cout << "Plain query result: " << std::endl;
    std::cout << std::setw(12) << "revenue" << std::endl;
    std::cout << std::setw(12) << plain_agg_res << std::endl;

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    return filtering_time + aggregation_time;
}

int main(int argc, char *argv[])
{ // std::chrono::system_clock::time_point start, end;
    // start = std::chrono::system_clock::now();
    // std::ifstream RQ1("../../engorgio_relational_q1_2.csv", std::ios::app);
    // std::vector<std::vector<std::string>> csvData_rq1;
    // std::string line;
    // while (std::getline(RQ1, line))
    // {
    //     std::vector<std::string> row;
    //     std::stringstream ss(line);
    //     std::string cell;
    //     while (std::getline(ss, cell, ','))
    //     {
    //         row.push_back(cell);
    //     }
    //     csvData_rq1.push_back(row);
    // }
    // RQ1.close();

    // for (int i = 2048; i < 16385; i *= 2)
    // {
    //     double rq1_time = Eval_E2E_q1(i);
    //     for (int j = 1; j < csvData_rq1.size(); j++)
    //     {
    //         if (std::stoi(csvData_rq1[j][0]) == i)
    //         {
    //             csvData_rq1[j][3] = std::to_string(rq1_time);
    //             break;
    //         }
    //     }
    // }
    // std::ofstream outFile_RQ1("../../engorgio_relational_q1_2.csv");
    // for (const auto &row : csvData_rq1)
    // {
    //     for (int i = 0; i < row.size(); i++)
    //     {
    //         outFile_RQ1 << row[i];
    //         if (i < row.size() - 1)
    //         {
    //             outFile_RQ1 << ",";
    //         }
    //     }
    //     outFile_RQ1 << std::endl;
    // }
    // outFile_RQ1.close();

    // std::ifstream RQ12("../../engorgio_relational_q12_2.csv", std::ios::app);
    // std::vector<std::vector<std::string>> csvData_rq12;
    // while (std::getline(RQ12, line))
    // {
    //     std::vector<std::string> row;
    //     std::stringstream ss(line);
    //     std::string cell;
    //     while (std::getline(ss, cell, ','))
    //     {
    //         row.push_back(cell);
    //     }
    //     csvData_rq12.push_back(row);
    // }
    // RQ12.close();
    // for (int i = 2048; i < 16385; i *= 2)
    // {

    //     double rq12_time = Eval_E2E_q12(i);
    //     for (int j = 1; j < csvData_rq12.size(); j++)
    //     {
    //         if (std::stoi(csvData_rq12[j][0]) == i)
    //         {
    //             csvData_rq12[j][3] = std::to_string(rq12_time);
    //             break;
    //         }
    //     }
    // }
    // std::ofstream outFile_RQ12("../../engorgio_relational_q12_2.csv");
    // for (const auto &row : csvData_rq12)
    // {
    //     for (int i = 0; i < row.size(); i++)
    //     {
    //         outFile_RQ12 << row[i];
    //         if (i < row.size() - 1)
    //         {
    //             outFile_RQ12 << ",";
    //         }
    //     }
    //     outFile_RQ12 << std::endl;
    // }
    // outFile_RQ12.close();
    // end = std::chrono::system_clock::now();
    // double total_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000;
    // std::cout << "RQ compute time:" << total_time << std::endl;
    // return 0;
    Eval_E2E_q6(65536);
}
