#ifndef clockType_hh
#define clockType_hh

class clockType
{
public:
	void setTime(int, int, int);
	void getTime(int&, int&, int&) const;
	void printTime() const;
	void incrementSeconds();
	void incrementMinutes();
	void incrementHours();
	bool equalTime(const clockType&) const;

private:
	int hr;
	int min;
	int sec;	
};

#endif // clockType_hh