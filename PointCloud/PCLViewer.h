#pragma once
class PCLViewer
{
public:
	PCLViewer();
	~PCLViewer();

	void ReadPcdFile(pcl::PointCloud<pcl::PointXYZRGBA>::Ptr cloud, const fs::path &fileName);		//	��ȡpcd�ļ��еõ�
private:
	void ReadBuffer(const fs::path &fileDir, char **buffer);
private:
	unsigned int red;
	unsigned int green;
	unsigned int blue;
};

