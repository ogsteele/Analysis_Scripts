%% Penn Peak Scale (python)
% 
% def peakscale():
%     """
%     Scale the selected traces in the currently active channel to their mean peak amplitude.
%     """
% 
%     # Measure baseline in selected traces
%     base=[]
%     for i in stf.get_selected_indices():
%         stf.set_trace(i)
%         base.append(stf.get_base())
% 
%     # Subtract baseline from selected traces
%     stf.subtract_base()
% 
%     # Measure peak amplitudes in baseline-subtracted traces
%     stf.select_all()
%     peak = []
%     for i in stf.get_selected_indices():
%         stf.set_trace(i)
%         peak.append(stf.get_peak())
% 
%     # Calculate scale factor to make peak equal to the mean peak amplitude
%     scale_factor = peak / np.mean(peak)
% 
%     # Scale the traces and apply offset equal to the mean baseline
%     scaled_traces = [stf.get_trace(i) / scale_factor[i] + np.mean(base) for i in stf.get_selected_indices()]
% 
%     # Close window of baseline-subtracted traces
%     stf.close_this()
% 
%     return stf.new_window_list(scaled_traces)
%     
%     
%     def subtract_base():
%     """
%     """
%     subtracted_traces = []
%     for i in range(stf.get_size_channel()):
%         stf.set_trace(i);
%         subtracted_traces.append(stf.get_trace() - stf.get_base())
%     stf.new_window_list(subtracted_traces)