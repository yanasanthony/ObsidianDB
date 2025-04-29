
[**Step1**] 获取所研究复杂疾病的样本以及测序的基因列表，确保各基因值已经经过TPM标准化以及进行相关预处理

---
[**Step2**] 将基因列表与STRING数据库中的PPI网络匹配，去除仅有单个节点而无任何邻接节点的基因

---

[**Step3**] 由于PPI网络映射得到的子网路是无方向性的，需要通过构建因果强度指数$w$来判断。定义计算时间$t$时一个病例样本的基因间因果强度：$w_{in}(g_{j}^{k} \rightarrow g^{k})=SD(g_{j}^{k})\cdot PCC(g_{j}^{k}\leftrightarrow g^{k})$
由PPI网络获得相互作为阈值作为互信息度：$I(g_{j}^{k};g^{k}):I(g_{k}^{k};g^{k})>\beta*max_{j}I(g_{j}^{k};g^{k})=cutoff$

---

[**Step4**] 从时间$t$时该病例样本的全局因果网络$N_{S}$中提取各局部因果网络$\{LN_{S}^{k}\}_{k=1}^{Q}$，局部因果网络又可分为局部入度因果网络$\{N_{S}^{k}\}_{k=1}^{Q}$与局部出度因果网络$\{L_{S}^{k}\}_{k=1}^{Q}$，定义由邻接基因指向中心基因为入度（in-degree），由中心基因指向邻接基因为出度（out-degree）。
- 中心基因为$g^{k}$的局部入度因果网络$N_{S}^{k}=\{g_{in,1}^{k},g_{in,2}^{k},...,g_{in,N}^{k}\}$，有与中心基因的因果强度$\{w_{in}(g_{in,1}^{k} \rightarrow g^{k}),w_{in}(g_{in,2}^{k} \rightarrow g^{k}),...,w_{in}(g_{in,N}^{k} \rightarrow g^{k})\}$
- 中心基因为$g^{k}$的局部出度因果网络$L_{S}^{k}=\{g_{out,1}^{k},g_{out,2}^{k},...,g_{out,L}^{k}\}$，有与中心基因的因果强度$\{w_{out}(g^{k} \rightarrow g_{out,1}^{k}),w_{out}(g^{k} \rightarrow g_{out,2}^{k}),...,w_{out}(g^{k} \rightarrow g_{out,L}^{k})\}$ 

---

[**Step5**] 计算各局部因果网络的***RF-SCNE***值
$$
H^{k}=\alpha H_{in}^{k}+(1-\alpha)H_{out}^{k},\alpha=\frac{\sum_{j=1}^{N}w_{in}(g_{j}^{k} \rightarrow g^{k})}{\sum_{j=1}^{N}w_{in}(g_{j}^{k} \rightarrow g^{k}) + \sum_{j=1}^{L}w_{out}(g^{k} \rightarrow g_{out}^{k})}
$$
$$
\left\{
\begin{array}{**lr**}
H_{in}^{k}=-\frac{\sum_{j=1}^{N}p_{in,j}log(p_{in,j})\sigma(FT(g^{k}))e^{w_{in,j}}}{N},FT(g^{k})=\left|\frac{E_{case}(g^{k})-\mu_{ref}(g^{k})}{\sigma_{ref}(g^{k})}\right|,\sigma(FT(g^{k}))=\frac{1}{1+e^{-FT(g^{k})}}, & \\
H_{out}^{k}=-\frac{\sum_{j=1}^{L}p_{out,j}log(p_{out,j}) \sigma( \overline{FT}(g^{k}))e^{w_{out,j}}}{L},\overline{FT}(g^{k})=\frac{\sum_{j=1}^{L}\left|\frac{E_{case}(g_{j}^{k})-\mu_{ref}(g_{j}^{k})}{\sigma_{ref}(g_{j}^{k})}\right|}{L},\sigma(\overline{FT}(g^{k}))=\frac{1}{1+e^{-\overline{FT}(g^{k})}}
\end{array}
\right.
$$
其中，加权概率分布如下：
$$
\left\{
\begin{array}{**lr**}
p_{in,j}=\frac{w_{in}(g_{j}^{k} \rightarrow g^{k})*I(g_{j}^{k};g^{k})}{\sum_{j=1}^{N}w_{in}(g_{j}^{k} \rightarrow g^{k})*I(g_{j}^{k};g^{k})}, & \\
p_{out,j}=\frac{w_{out}(g^{k} \rightarrow g_{j}^{k})*I(g^{k};g_{j}^{k})}{\sum_{j=1}^{N}w_{out}(g^{k} \rightarrow g_{j}^{k})*I(g^{k};g_{j}^{k})}
\end{array}
\right.
$$

---

[**Step6**] 计算全局因果网络的***RF-SCNE***值
$$
H_{t}=\sum_{k=1}^{Q}H^{k},Q=top 5\% * N(局部因果网络的RF-SCNE值序后的前5\%中心基因)
$$

---

[**Step7**] 使用配对样本$t$检验（如样本量大，采用正态检验）以判断时间点$t$较时间点$t-1$的$Q$个基因的表达量是否显著增大。
- 假设：
$$
\left\{
\begin{array}{**lr**}
H_{0}:时间点t与时间点t-1的Q个基因的表达量并未显著差别, & \\
H_{1}:时间点t较时间点t-1的Q个基因的表达量显著增大, & \\
H_{2}:时间点t较时间点t-1的Q个基因的表达量显著减小.
\end{array}
\right.
$$
- 计算统计量：用时间点$t$时的表达量$\{E_{t}^{k}\}_{k=1}^{Q}$与时间点$t-1$时的表达量$\{E_{t-1}^{k}\}_{k=1}^{Q}$做差，得到$d$值，计算$d$值的均数$\overline{d}$与标准差$S_{d}$，$t$统计量如下：
$$
t=\frac{\overline{d}}{S_{d}/ \sqrt{Q}}
$$
- 计算$P$值：在本案例中使用双侧检验，$t_{\frac{\alpha}{2}}(Q-1)$或$z_{\frac{\alpha}{2}}$为正临界值，结论如下：
$$
\left\{
\begin{array}{**lr**}
如果t>t_{\frac{\alpha}{2}}(Q-1),则拒绝H_{0},接受H_{1},认为显著增大, & \\
如果-t_{\frac{\alpha}{2}}(Q-1) \leq t \leq t_{\frac{\alpha}{2}}(Q-1),则接受H_{0},拒绝H_{1},H_{2},认为显著增大, & \\
如果t<-t_{\frac{\alpha}{2}}(Q-1),则拒绝H_{0},接受H_{2},认为显著减小
\end{array}
\right.
$$
- 做出结论