Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())

getwd()
setwd('D:/GEO/video20230606/KEGG pathway')



library(ggplot2)
library(RColorBrewer)

dt = read.delim(file.choose(),row.names = 1)

colnames(dt)

dt$Term <- factor(dt$Term,
                  levels = rev(dt$Term))

mytheme <- theme(
  axis.title = element_text(size = 13),
  axis.text = element_text(size = 11),
  plot.title = element_text(size = 14,
                            hjust = 0.5,
                            face = "bold"),
  legend.title = element_text(size = 13),
  legend.text = element_text(size = 11),
  plot.margin = margin(t = 5.5,
                       r = 10,
                       l = 5.5,
                       b = 5.5)
)

p <- ggplot(data = dt, aes(x = Count, y = Term, fill = -log10(PValue))) +
  scale_fill_distiller(palette = "RdPu",direction = 1) + #������ɫ
  geom_bar(stat = "identity", width = 0.8) + #��������ͼ
  labs(x = "Number of Gene", y = "", title = "KEGG enrichment barplot") + #�޸�/���Ӹ�����
  theme_bw() + mytheme #�������
p



mytheme2 <- mytheme + theme(axis.text.y = element_blank())
p1 <- ggplot(data = dt, aes(x = Fold.Enrichment, y = Term, fill = -log10(PValue))) +
  scale_fill_gradient(low = "tomato", high = "brown4") + #����ɫtomato��brown4����
  geom_bar(stat = "identity", width = 0.8, alpha = 0.9) + #bar�Ŀ��Ⱥ�͸����width /alpha. ��������ͼ��ࣺ���Ϻ��� position = position_dodge(width = 0.5)
  labs(x = "RichFactor", y = "", title = "GO MF enrichment term") +
  geom_text(aes(x = 0.03, #����ֵ���������ı���ǩ��ʼλ��
                label = Term),
            hjust = 0, size = 4)+ #hjust = 0,�����,size���������С
  theme_classic() + mytheme2

p1
