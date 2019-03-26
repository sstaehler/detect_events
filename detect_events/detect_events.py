#!/usr/bin/env python
"""

"""
__author__ = "Simon Staehler"

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from obspy.signal.util import next_pow_2
from obspy import read
import instaseis
from obspy import UTCDateTime as utct

def pickify(fig, ax, dist, M, t0):
    keymap = {'w': 0, 'e': 1}
    def onkey(event):
        global events
        if event.key in keymap.keys():
            events.append([dist, M, keymap[event.key]])
        elif event.key == 'r':
            print('removed event from list')
            del events[-1]
            events.append([dist, M, 0])
            plt.close()

        ax.axvline(x=t0 + dist * (66 / 2.5), color='lightgreen', lw=2)
        ax.axvline(x=t0 + dist * (66 / 3.5), color='lightgreen', lw=2)

        fig.suptitle(('Event: Magnitude %3.1f in %3d degree distance. '
                     ' Press "q" to save, "r" to remove (in case of wrong '
                      'pick)')
                     % (M, dist))
        fig.canvas.draw()


    fig.canvas.mpl_connect('key_press_event', onkey)

    plt.show()

    return


def plot_spec(st_HF, st_LF, winlen_sec_LF=200, winlen_sec_HF=10., overlap=0.5):
    fig = plt.figure(figsize=(16, 8))
    # [left bottom width height]
    ax_seis_LF = fig.add_axes([0.06, 0.78, 0.8, 0.2], label='seismogram LF')
    ax_seis_HF = ax_seis_LF.twinx()
    ax_spec_LF = fig.add_axes([0.06, 0.13, 0.8, 0.4], sharex=ax_seis_LF,
                              label='spectrogram LF')
    ax_spec_HF = fig.add_axes([0.06, 0.53, 0.8, 0.25], sharex=ax_seis_LF,
                              label='spectrogram HF')

    winlen_LF = winlen_sec_LF * st_LF[0].stats.sampling_rate
    winlen_HF = winlen_sec_HF * st_HF[0].stats.sampling_rate

    freq_HF = st_HF[0].stats.sampling_rate
    fmax_LF = 1.0
    fmax_HF = freq_HF / 2.
    fmin_LF = 2. / (winlen_sec_LF)
    fmin_HF = 1.0
    # Colorbar axis
    ax_cb = fig.add_axes([0.92, 0.8, 0.03, 0.15], label='colorbar')

    for tr in st_LF:
        t0 = float(tr.stats.starttime)
        ax_seis_LF.plot((tr.times() + t0), tr.data * 1e9,
                        'navy', lw=0.3)
    for tr in st_HF:
        t0 = float(tr.stats.starttime)
        ax_seis_HF.plot((tr.times() + t0), tr.data * 1e9,
                        'darkred', lw=0.3)

    for st, ax_spec, flim, winlen in zip([st_LF, st_HF],
                                                 [ax_spec_LF, ax_spec_HF],
                                                 [(fmin_LF, fmax_LF),
                                                  (fmin_HF, fmax_HF)],
                                                 [winlen_LF, winlen_HF]):
        for tr in st:
            p, f, t = mlab.specgram(tr.data, NFFT=winlen,
                                    Fs=tr.stats.sampling_rate,
                                    noverlap=winlen * overlap,
                                    pad_to=next_pow_2(winlen) * 4)
            t += float(tr.stats.starttime)

            bol = np.array((f > flim[0], f < flim[1])).all(axis=0)

            vmin = 1e-20 #np.percentile(p[bol, :], q=1, axis=None)
            vmax = 1e-15 #np.percentile(p[bol, :], q=90, axis=None)

            ax_spec.pcolormesh(t, f[bol], 10 * np.log10(p[bol, :]),
                               vmin=10 * np.log10(vmin),
                               vmax=10 * np.log10(vmax),
                               cmap='plasma')

    ax_seis_LF.set_ylim(st_LF[0].std() * 1e10 * np.asarray([-2., 1.]))
    ax_seis_HF.set_ylim(st_HF[0].std() * 1e10 * np.asarray([-1., 2.]))
    ax_seis_LF.set_ylabel('%s (<1Hz) [nm/s]' % st_LF[0].stats.channel,
                          color='navy')
    ax_seis_HF.set_ylabel('%s (>1Hz) [nm/s]' % st_HF[0].stats.channel,
                          color='darkred')
    ax_seis_LF.tick_params('y', colors='navy')
    ax_seis_HF.tick_params('y', colors='darkred')
    ax_spec_HF.set_ylim(fmin_HF, fmax_HF)
    ax_spec_HF.set_ylabel('frequency / Hz', fontsize=12)

    ax_spec_LF.set_yscale('log')
    ax_spec_LF.set_ylabel('period / seconds', fontsize=12)

    # This needs to be done for both axes, otherwise the PSD and the Spec axis
    # plot weird yticks on top of each other
    for ax in [ax_spec_LF]:
        ax.set_ylim(fmin_LF, fmax_LF)
        tickvals = np.asarray([1./2, 1./5, 1./10, 1./20, 1./50, 1./100, 1./200.])
        ax.set_yticks(tickvals[tickvals > fmin_LF])
        ax.set_yticklabels(['2', '5', '10', '20', '50', '100', '200'])
        ax.set_yticklabels([], minor=True)

    # make unnecessary labels disappear
    for ax in [ax_spec_HF, ax_seis_LF]:
        plt.setp(ax.get_xticklabels(), visible=False)

    # Axis with colorbar
    mappable = ax_spec_HF.collections[0]
    cb = plt.colorbar(mappable=mappable, cax=ax_cb)
    ax_cb.set_ylabel('PSD (m/s)/Hz')
    fig.suptitle('Press "e" for event detected, "w" for none detected')


    return fig, ax_spec_LF

def create_event(M_min=2.0, M_max=4.0):
    M = np.random.rand(1) * (M_max - M_min) + M_min
    dist = 2 + (np.random.rand(1)**2) * 178.
    return M, dist


def autoreject(M, dist):
    reject = False
    if M<2.5 and dist > 30:
        reject = True
    return reject

def autoaccept(M, dist):
    reject = False
    if M>3.0 and dist < 20:
        reject = True
    return reject

def select_and_add(st, db_HF, db_LF, M, dist):
    winlen_hours = 4
    tstart_total = float(st[0].stats.starttime) + 3600
    tend_total = float(st[0].stats.endtime) - 3600 * (winlen_hours + 1)

    tstart_win = tstart_total + (tend_total - tstart_total) * np.random.rand(1)
    tend_win = tstart_win + 3600 * winlen_hours

    t0 = tstart_win + 600 + (tend_win - 4200 - tstart_win) * np.random.rand(1)

    st_LF = st.slice(starttime=utct(tstart_win - 1800),
                     endtime=utct(tend_win + 1800))
    st_HF = st.slice(starttime=utct(tstart_win - 1800),
                     endtime=utct(tend_win + 1800))

    st_HF.filter('highpass', freq=0.9)
    st_LF.filter('bandpass', freqmin=1./60, freqmax=1.2)

    st_LF.trim(starttime=utct(tstart_win),
               endtime=utct(tend_win))
    st_HF.trim(starttime=utct(tstart_win),
               endtime=utct(tend_win))

    src = instaseis.Source(latitude=90.0-dist,
                           longitude=0.0,
                           depth_in_m=1e3 + np.random.rand(1) * 2.9e4,
                           m_rr=instaseis.source.magnitude2moment(M)[0] *
                                np.sqrt(2.),
                           origin_time=utct(t0)
                           )
    rec = instaseis.Receiver(latitude=90.0, longitude=0.0,
                             network='XB', station='ELYSE', location='58')

    for db, flims in zip([db_HF, db_LF], [(1./10, 4.), (1./100., 1./10.)]):
        st_inst = db.get_seismograms(source=src, receiver=rec,
                                     kind='velocity',
                                     components='Z', dt=0.1)
        st_inst.trim(endtime=st_HF[0].stats.endtime-2)
        st_inst[0].stats.channel = 'BZC'
        #st_inst.filter('lowpass', freq=fmin)
        st_inst.filter('highpass', freq=flims[0])
        st_inst.filter('lowpass', freq=flims[1])
        i0 = int((st_inst[0].stats.starttime - st_LF[0].stats.starttime) * 10)
        for s in (st_LF, st_HF):
            s[0].data[i0:i0+len(st_inst[0].data)] += st_inst[0].data

    return st_LF, st_HF, t0


def play(path_LF, path_HF, nevent=40):
    global events
    events = []
    st = read('VBB_3days_VEL.mseed')

    db_HF = instaseis.open_db(path_HF)
    db_LF = instaseis.open_db(path_LF)
    for i in range(0, nevent):
        M, dist = create_event(M_min=2.0, M_max=4.0)
        if autoreject(M, dist):
            events.append([dist, M, 0])
            string = '%3d/%3d: Below threshold: Event in %d deg with M%3.1f'
        elif autoaccept(M, dist):
            events.append([dist, M, 1])
            string = '%3d/%3d: Above threshold: Event in %d deg with M%3.1f'
        else:
            st_LF, st_HF, t0 = select_and_add(st, db_HF, db_LF, M, dist)
            fig, ax = plot_spec(st_HF=st_HF, st_LF=st_LF, winlen_sec_LF=50)
            picks = pickify(fig, ax, M=M, dist=dist, t0=t0)
            if events[-1][-1] == 1:
                string = '%3d/%3d: Found event in %d deg with M%3.1f'
            else:
                string = '%3d/%3d: Missed event in %d deg with M%3.1f'
        print(string % (i, nevent, dist, M))
    with open('result.txt', 'a') as f:
        for event in events:
            plt.plot(event[0], event[1], 'o', c='C%d' % event[2])
            f.write('%5.1f,  %4.2f, %d \n' % (event[0], event[1], event[2]))
    plt.show()


if __name__ == '__main__':
    nevent = 40

    # The good ones
    path_LF = 'http://instaseis.ethz.ch/blindtest_5s/EH45TcoldCrust1b_5s'
    path_HF = 'http://instaseis.ethz.ch/blindtest_1s/EH45TcoldCrust1b_Q100_1s'
    play(path_LF, path_HF, nevent)