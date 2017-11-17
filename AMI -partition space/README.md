
Adjusted Mutual Information between multi-scale/multi-resolution partition solution found by stability analysis
===============================================================================================================

Distance function that provided the best results is based on the mutual information ($MI$) between stress-specific partition solutions at different resolutions. Yet, it is corrected for chance under an hyper-geometric model (Eq.~\ref{eq:Eq_AMI}) by the expected value for $MI$ between two partitions, $\mathbb{E}[MI(\Gamma_{t}^{Si},\Gamma_{t'}^{Sj})]$, and the entropy of each partition, $H(\Gamma_{t}^{Si})$.

\begin{equation}
\begin{split}
&AMI(\Gamma_{t}^{Si},\Gamma_{t'}^{Sj})=\\
&\frac{MI(\Gamma_{t}^{Si},\Gamma_{t'}^{Sj})-\mathbb{E}[MI(\Gamma_{t}^{Si},\Gamma_{t'}^{Sj})]}{max(H(\Gamma_{t}^{Si}),H(\Gamma_{t'}^{Sj}))-\mathbb{E}[MI(\Gamma_{t}^{Si},\Gamma_{t'}^{Sj})]}\label{eq:Eq_AMI}
\end{split}
\end{equation} 

References on AMI:

[1] 'A Novel Approach for Automatic Number of Clusters Detection based on Consensus Clustering', N.X. Vinh, and Epps, J., in Procs. IEEE Int. Conf. on Bioinformatics and Bioengineering (Taipei, Taiwan), 2009.

[2] 'Information Theoretic Measures for Clusterings Comparison: Is a Correction for Chance Necessary?', N.X. Vinh, Epps, J. and Bailey, J., in Procs. the 26th International Conference on Machine Learning (ICML'09)

[3] 'Information Theoretic Measures for Clusterings Comparison: Variants, Properties, Normalization and Correction for Chance', N.X. Vinh, Epps, J. and Bailey, J., Journal of Machine Learning Research, 11(Oct), pages 2837-2854, 2010
