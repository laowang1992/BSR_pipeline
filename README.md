# BSR_pipeline
a pipeline of BSR-seq, without parents information
# 原理
BSR-seq属于BSA分析方法中的一种，只不过是使用转录组测序的数据进行分析。在这里只有两个混池的转录组测序数据而缺少亲本数据，因此使用欧几里得距离（欧氏距离，Euclidean Distance，ED）进行分析，即：![](http://latex.codecogs.com/gif.latex?ED=\sqrt{(A_{mut}-A_{wt})^2+(C_{mut}-C_{wt})^2+(G_{mut}-G_{wt})^2+(G_{mut}-G_{wt})^2+(T_{mut}-T_{wt})^2})。同时，因为没有亲本信息，在SNP位点过滤也需要注意，首先去掉低质量位点，然后去掉两池纯和且相同的位点和存在缺失值的位点。在最后画图时可以考虑使用ED<sup>4</sup>以放大信号。
# 结果
## 查看SNP分布
<img src="https://user-images.githubusercontent.com/35584208/125470830-31c81962-9fc1-4b59-802b-4c5921ae27ea.png" width = "400" alt="SNP_distribution_histogram" /><br>
## ED plot
![BSR_ED](https://user-images.githubusercontent.com/35584208/125471595-6d7877a1-2b07-4177-b60c-5b745e79d12f.png)
## ED<sup>4</sup> plot
![BSR_ED4](https://user-images.githubusercontent.com/35584208/125471845-246f4bc9-5f22-4508-a198-6eba1c24789e.png)


