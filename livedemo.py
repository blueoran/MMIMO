import os
import numpy as np
import trace.Decoder
from tkinter import *
from tkinter import messagebox
from subprocess import Popen, PIPE

mod_symbol=[]

def gen_symbols(mod_order):
    symbols = []
    num = 0
    if mod_order >= 4:
        for i in (-1, 1):
            for j in (-1, 1):
                symbols.append(complex(i, j))
    if mod_order >= 16:
        for j in (-1, 1):
            symbols.append(complex(-3, j))
        for j in (-1, 1):
            symbols.append(complex(3, j))
        for i in (-1, 1):
            symbols.append(complex(i, -3))
        for i in (-1, 1):
            symbols.append(complex(i, 3))
        for i in (-3, 3):
            for j in (-3, 3):
                symbols.append(complex(i, j))
    if mod_order >= 32:
        for j in range(-3, 4, 2):
            symbols.append(complex(-5, j))
        for j in range(-3, 4, 2):
            symbols.append(complex(5, j))
        for i in range(-3, 4, 2):
            symbols.append(complex(i, -5))
        for i in range(-3, 4, 2):
            symbols.append(complex(i, 5))
    if mod_order >= 64:
        for i in (-5, 5):
            for j in (-5, 5):
                symbols.append(complex(i, j))
        for j in range(-5, 6, 2):
            symbols.append(complex(-7, j))
        for j in range(-5, 6, 2):
            symbols.append(complex(7, j))
        for i in range(-5, 6, 2):
            symbols.append(complex(i, -7))
        for i in range(-5, 6, 2):
            symbols.append(complex(i, 7))
        for i in (-7, 7):
            for j in (-7, 7):
                symbols.append(complex(i, j))
    return symbols

def random_gen():
    # messagebox.showinfo('message', f'give flower{int(random_gen_entry_num.get())}{ofdm_num[ofdm_symbol.get()]}')  # 提示框
    send_data.delete(1.0,END)
    decoder_type_str=decoder_type.get()
    is_verbose_yesno=is_verbose.get()
    ofdm_symbol_number=ofdm_symbol.get()
    sim_times_number=sim_times.get()
    num_sender_number=num_sender.get()
    num_receiver_number=num_receiver.get()
    add_noise_yesno=add_noise.get()
    log10_SNR_number=log10_SNR.get()
    # print(ofdm_symbol_num)
    if num_sender_number>num_receiver_number:
        send_data.delete(1.0,END)
        send_data.insert(END,"Receiver's Antennas Number should not less than Sender's Antennas Number!!")
        return
    if_verbose='OPT+=-DVERBOSE' if is_verbose_yesno else''
    os.system(f'cd {os.path.dirname(__file__)}; make clean; make {decoder_type_str} {if_verbose}')
    # os.system(f'{os.path.dirname(__file__)}/mysim.exe {ofdm_symbol_number} {sim_times_number} {num_sender_number} {num_receiver_number} {add_noise_yesno} {log10_SNR_number}')
    p = Popen(f'{os.path.dirname(__file__)}/mysim.exe {ofdm_symbol_number} {sim_times_number} {num_sender_number} {num_receiver_number} {add_noise_yesno} {log10_SNR_number}'.split(' '), stdout=PIPE, stderr=PIPE, stdin=PIPE)
    output = p.stdout.read()
    send_data.insert(END,bytes.decode(output))


def simulation():
    pass

ofdm=["QPSK","16-QAM","32-QAM","64-QAM"]
ofdm_num=[4,16,32,64]
ofdm_symbols={ofdm_num[i]:gen_symbols(ofdm_num[i]) for i in range(4)}

'''
" <module num> <sim times> <sender num> <receiver num> <add "
                "noise> <log10-SNR>" 
'''
root=Tk()
root.title("LiveDemo for our MMIMO Decoder")

decoder_type=StringVar()
label_decoder_type=Label(root,text='Decoder Type:')
label_decoder_type.pack()
decoder_types=['Zero Forcing','Sphere Decoding','Sphere Decoding with Radius Initializing Optimization','K-Best Sphere Decoding']
decoder_cmds=['OPT+=-DZF','OPT+=-DSP','OPT+=-DSP OPT+=-DSP_RADIUS_OPT','OPT+=-DKSD']
decoder_radio_button=[Radiobutton(root,text=decoder_types[i],variable=decoder_type,value=decoder_cmds[i]) for i in range(4)]
for i in range(4):
    decoder_radio_button[i].pack()

is_verbose=IntVar()
is_verbose.set(0)
check_is_verbose=Checkbutton(root, text="Verbose",variable=is_verbose,onvalue=1, offvalue=0)
check_is_verbose.pack()

# random_gen_entry_num=StringVar()
label_ofdm_symbol=Label(root,text='OFDM Symbol:')
label_ofdm_symbol.pack()
ofdm_symbol=IntVar()
ofdm_symbol.set(32)
ofdm_radio_button=[Radiobutton(root,text=ofdm[i],variable=ofdm_symbol,value=ofdm_num[i]) for i in range(4)]
for i in range(4):
    ofdm_radio_button[i].pack()
    
label_sim_times=Label(root,text='Simulation Times:')
label_sim_times.pack()
sim_times=IntVar()
sim_times.set(4)
num_sym_times=Entry(root,textvariable=sim_times)
num_sym_times.pack()

    
label_sender=Label(root,text='Sender Antennas Number:')
label_sender.pack()
num_sender=IntVar()
num_sender.set(4)
num_sender_entry=Entry(root,textvariable=num_sender)
num_sender_entry.pack()

label_receiver=Label(root,text='Receiver Antennas Number:')
label_receiver.pack()
num_receiver=IntVar()
num_receiver.set(8)
num_receiver_entry=Entry(root,textvariable=num_receiver)
num_receiver_entry.pack()

add_noise=IntVar()
add_noise.set(0)
check_add_noise=Checkbutton(root, text="Add Noise",variable=add_noise,onvalue=1, offvalue=0)
check_add_noise.pack()

log10_SNR=DoubleVar()
log10_SNR.set(0)
label_log10_SNR=Label(root,text='log_10(SNR):')
label_log10_SNR.pack()
num_log10_SNR=Entry(root,textvariable=log10_SNR)
num_log10_SNR.pack()




random_gen_button = Button(root,command=random_gen,text='Start Simulation')
random_gen_button.pack()

send_data=Text(root)
send_data.pack()


# simulation_button=Button(root,command=random_gen,text='Generate Random Symbols')

    

# print(gen_symbols(16))
# bt.bind('<Button-1>', random_gen)  # 绑定点击事件
root.mainloop()
