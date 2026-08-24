import functools

def verify_rns_comparison():


    rns_moduli = [3, 5, 7] 
    

    M = functools.reduce(lambda x, y: x * y, rns_moduli)
    
    print(f"--- RNS 设置 ---")
    print(f"模数基: {rns_moduli}")
    print(f"动态范围 M: {M} (奇偶性: {'奇数' if M % 2 != 0 else '偶数'})")
    print("-" * 30)


    if M % 2 == 0:
        print("错误：M 必须是奇数，请修改模数基。")
        return



    correct_count = 0
    total_count = 0
    


    for A in range(M):
        for B in range(M):
            total_count += 1
            




            D_rns = (A - B) % M
            

            lsb_a = A % 2
            lsb_b = B % 2
            lsb_d = D_rns % 2
            





            



            flip = lsb_a ^ lsb_b ^ lsb_d
            

            if flip == 1:
                inferred_result = "A < B"
            else:
                inferred_result = "A >= B"
                

            if A < B:
                actual_result = "A < B"
            else:
                actual_result = "A >= B"
            

            if inferred_result != actual_result:
                print(f"验证失败! A={A}, B={B}, D_rns={D_rns}")
                print(f"LSB(A)={lsb_a}, LSB(B)={lsb_b}, LSB(D)={lsb_d}, Flip={flip}")
                return
            else:
                correct_count += 1

    print(f"验证完成。")
    print(f"测试用例总数: {total_count}")
    print(f"通过用例总数: {correct_count}")
    
    if correct_count == total_count:
        print("\n结论: 逻辑验证成功！")
        print("公式: Is_Less_Than = LSB(A) XOR LSB(B) XOR LSB(D_rns)")
    else:
        print("\n结论: 逻辑存在漏洞。")


def demo_single_case(A, B, M):
    print(f"\n[单步演示] A={A}, B={B}, M={M}")
    D_rns = (A - B) % M
    
    lsb_a = A % 2
    lsb_b = B % 2
    lsb_d = D_rns % 2
    
    expected_arithmetic_lsb = lsb_a ^ lsb_b
    
    print(f"1. 计算 RNS 差值: ({A} - {B}) mod {M} = {D_rns}")
    print(f"2. 提取 LSB: LSB(A)={lsb_a}, LSB(B)={lsb_b}")
    print(f"3. 理论算术差 LSB (A-B): {lsb_a} ^ {lsb_b} = {expected_arithmetic_lsb}")
    print(f"4. 实际 RNS 差 LSB (D): {lsb_d}")
    
    flip = expected_arithmetic_lsb ^ lsb_d
    print(f"5. 比较差异 (XOR): {expected_arithmetic_lsb} ^ {lsb_d} = {flip}")
    
    if flip == 1:
        print("-> 结果: 检测到 LSB 翻转，判定 A < B")
    else:
        print("-> 结果: LSB 一致，判定 A >= B")

if __name__ == "__main__":
    verify_rns_comparison()
    

    # M = 105

    demo_single_case(10, 3, 105)

    demo_single_case(3, 10, 105)