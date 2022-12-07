# MMIMO

## 接口

```cpp
./decoder.cpp
complex<double> **Decoder(
                        int mod_order, 
                        int num_sender, 
                        int num_receiver, 
                        int num_ofdm_sym, 
                        complex<double> **H, 
                        complex<double> **Y, 
                        complex<double> *w = nullptr
                        )
```

$$
Y=Hx+w
$$
|param|type|description|
|:---:|:---:|:---:|
|mod_order|int|调制阶数|
|num_sender|int|发送天线数|
|num_receiver|int|接收天线数|
|num_ofdm_sym|int|OFDM并发符号数|
|H|[num_receiver*num_sender]|信道矩阵|
|Y|[num_ofdm_sym*num_receiver]|接收信号|
|w|[num_receiver]|噪音|

## 测试

### Simulation

- environment:
  - matlab
  - g++

```bash
cd sim/RenewLab
matlab
rl_ofdm_mimo
(显示: MIMO algorithm(ZF, Conj, C): ) 输入: C
```

### Traces

- environment:
  - numpy
  - scipy
  - scikit-commpy
  - scikit-learn
  - h5py
  - matplotlib

```bash
cd trace/
python test.py
```

## TODO

- rl_ofdm_mmimo_sim.m
- Agora
- armadillo(等矩阵分析lib)
- better simulator for python
