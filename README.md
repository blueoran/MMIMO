# MMIMO

## �ӿ�

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
|mod_order|int|���ƽ���|
|num_sender|int|����������|
|num_receiver|int|����������|
|num_ofdm_sym|int|OFDM����������|
|H|[num_receiver*num_sender]|�ŵ�����|
|Y|[num_ofdm_sym*num_receiver]|�����ź�|
|w|[num_receiver]|����|

## ����

### Simulation

- environment:
  - matlab
  - g++

```bash
cd sim/RenewLab
matlab
rl_ofdm_mimo
(��ʾ: MIMO algorithm(ZF, Conj, C): ) ����: C
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

### Mysim

 - environment:
   - g++
   - eigen3
 
��`MMIMO/`��ʹ��`make`��Ĭ��ʹ��zero forcing��

ͨ��`make OPT+=[OPTION]`��ӱ���ѡ�

ͨ��`OPT+=-DSP`ʹ��sphere decoding��

ͨ��`OPT+=-DFSD`ʹ��fixed sphere decoding��FSD����

ͨ��`OPT+=-DSP_RADIUS_OPT`���ø�˹�����ķ������ȷ�������뾶�ĳ�ֵ����ʽ��$r^2 = \alpha n_{receiver} \sigma^2$��sphere decoding��FSD����ʹ�ô�ѡ�

ͨ��`OPT+=-DVERBOSE`�����ϸ�������

## TODO

- rl_ofdm_mmimo_sim.m
- Agora
- armadillo(�Ⱦ������lib)
- better simulator for python
